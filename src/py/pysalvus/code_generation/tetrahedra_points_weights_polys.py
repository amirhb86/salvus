from sympy import *
from sympy.utilities.codegen import CCodeGen
import numpy as np
import pylab as pl
import numpy.linalg
import sys
from IPython.core.debugger import Tracer

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

from BuildNodesTetrahedraP3 import BuildNodesTetrahedraP3

def p3QuadratureWeights():
    
    # reference weights from MJS Chin-Joe-Kong et al (1999)
    # weights are given based on reference triangle with V=1/6, which need to be normalized
    w0 = 0.2143608668049743e-03*6 # corners
    wa = 0.8268179517797114e-03*6 # edges
    wb1 = 0.1840177904191860e-02*6 # faces1
    wb2 = 0.1831324329245650e-02*6 # faces2
    wg = 0.7542468904648131e-02*6 # interior1
    wd = 0.1360991755970793e-01*6 # interior2

    # pts given in order:
    # {interior, faces, edges, vertices}
    
    quadrature_weights =  [wg, wg, wg, wg,
                           wd, wd, wd, wd, wd, wd,
                           wb1, wb1, wb1,
                           wb2, wb2, wb2,
                           wb1, wb1, wb1,
                           wb2, wb2, wb2,
                           wb1, wb1,wb1,
                           wb2, wb2, wb2,
                           wb1, wb1, wb1,
                           wb2, wb2, wb2,
                           wa, wa, wa, wa,
                           wa, wa, wa, wa,
                           wa, wa, wa, wa,
                           w0, w0, w0, w0]

    return quadrature_weights
        
def p3ReferenceNodes():
    
    # P3 reference tet in PETSc ordering (interior,face,edges,vertices)
    (r,s,t) = BuildNodesTetrahedraP3(plot=False)
    
    return (r,s,t)

def p3ElementPolynomials(rstn,r,s,t):
    # faces
    b_1 = r*s*(1-r-s-t)
    b_2 = r*t*(1-r-s-t)
    b_3 = t*s*(1-r-s-t)
    b_4 = r*s*t
    P2_b_1 = b_1 * np.asarray([r,s,r*s,r**2,s**2])
    P2_b_2 = b_2 * np.asarray([r,t,r*t,r**2,t**2])
    P2_b_3 = b_3 * np.asarray([t,s,t*s,t**2,s**2])
    rr = (1-r-t)
    ss = (1-s-t)
    P2_b_4 = b_4 * np.asarray([rr,ss,rr*ss,rr**2,ss**2])
    
    # interior
    btilde = r*s*t*(1-r-s-t)
    P1t = np.asarray([1,r,s,t,r*s,r*t,s*t,r**2,s**2,t**2])
    P1t_bt = btilde * P1t

    # edges
    P3 = np.asarray([1+0*r,
                     r,r**2,r**3,s,s**2,s**3,t,t**2,t**3,
                     r*s,r**2*s,r*s**2,
                     r*t,r**2*t,r*t**2,
                     s*t,s**2*t,s*t**2,
                     r*s*t])
    
    # full polyspace
    P3all = np.concatenate((P3,P1t_bt,P2_b_1,P2_b_2,P2_b_3,P2_b_4))
    
    Np=len(P3all)
    V = np.zeros((Np,Np))
    (rn,sn,tn) = rstn
    for (i,(ri,si,ti)) in enumerate(zip(rn,sn,tn)):        
        V[i,:] =  np.asarray([P3all_j.subs({r:ri,s:si,t:ti}) for P3all_j in P3all])

    condV = np.linalg.cond(V)
    print("Condition number of V = %e" % (condV))

    # matrix can't be singular
    assert condV < 1/sys.float_info.epsilon, "Matrix is likely singular, no Lagrange polynomial possible!"

    # build polynomial coefficients
    P_a = []
    for i in range(0,Np):
        b = np.zeros((Np,1))
        b[i] = 1
        a = numpy.linalg.solve(V,b)
        P_a.append(a)

    phi = []
    for (i,p_ai) in enumerate(P_a):
        phi_i = np.dot(P3all, p_ai)[0]
        phi.append(phi_i)

    return phi

def p3ElementDerivatives(phi,r,s,t):
    deriv_r = []
    deriv_s = []
    deriv_t = []
    for phi_i in phi:
        dphi_dr = diff(phi_i,r)
        dphi_ds = diff(phi_i,s)
        dphi_dt = diff(phi_i,t)
        deriv_r.append(dphi_dr)
        deriv_s.append(dphi_ds)
        deriv_t.append(dphi_dt)
    return (deriv_r,deriv_s,deriv_t)

# def genK(phi,wn,r,s,t,rstn):
#     Np = len(wn)
#     Kwn_rr = np.zeros((Np,Np))
#     Kwn_rs = np.zeros((Np,Np))
#     Kwn_sr = np.zeros((Np,Np))
#     Kwn_ss = np.zeros((Np,Np))
#     if os.path.isfile("p3_triangle_Krr_Krs_Ksr_Kss.npz"):
#         print("Using precomputed Krr,etc")
#         Kall = np.load("Krr_Krs_Ksr_Kss.npz")
#         Kwn_rr = Kall["Kwn_rr"]
#         Kwn_rs = Kall["Kwn_rs"]
#         Kwn_sr = Kall["Kwn_sr"]
#         Kwn_ss = Kall["Kwn_ss"]
#     else:
#         print("Computing Krr,etc and storing them")
#         for i in range(0,Np):
#             print("i={}".format(i))
#             for j in range(0,Np):
#                 # exact version
#                 # phiI_x_phiJ_rr = diff(phi[i],r)*diff(phi[j],r)
#                 # phiI_x_phiJ_ss = diff(phi[i],s)*diff(phi[j],s)
#                 # Kexact[i,j] = integrate(phiI_x_phiJ_rr,(r,-1,-s),(s,-1,1)) + integrate(phiI_x_phiJ_ss,(r,-1,-s),(s,-1,1))
#                 # phiI_x_phiJ_rr = diff(phi[i],r)*diff(phi[j],r)
#                 # phiI_x_phiJ_ss = diff(phi[i],s)*diff(phi[j],s)
#                 # phiI_x_phiJ_rs = diff(phi[i],r)*diff(phi[j],s)
#                 # phiI_x_phiJ_sr = diff(phi[i],s)*diff(phi[j],r)
#                 # Kwn_rr[i,j] = integrate(phiI_x_phiJ_rr,(r,-1,-s),(s,-1,1))
#                 # Kwn_rs[i,j] = integrate(phiI_x_phiJ_rs,(r,-1,-s),(s,-1,1))
#                 # Kwn_sr[i,j] = integrate(phiI_x_phiJ_sr,(r,-1,-s),(s,-1,1))
#                 # Kwn_ss[i,j] = integrate(phiI_x_phiJ_ss,(r,-1,-s),(s,-1,1))

#                 # version using quadrature
#                 phiI_r = lambdify((r,s),diff(phi[i],r))
#                 phiI_s = lambdify((r,s),diff(phi[i],s))
#                 phiJ_r = lambdify((r,s),diff(phi[j],r))
#                 phiJ_s = lambdify((r,s),diff(phi[j],s))
#                 Kwn_rr[i,j] = np.asarray([phiI_r(rn,sn)*phiJ_r(rn,sn) for (rn,sn) in rsn]).dot(p3QuadratureWeights())
#                 Kwn_rs[i,j] = np.asarray([phiI_r(rn,sn)*phiJ_s(rn,sn) for (rn,sn) in rsn]).dot(p3QuadratureWeights())
#                 Kwn_sr[i,j] = np.asarray([phiI_s(rn,sn)*phiJ_r(rn,sn) for (rn,sn) in rsn]).dot(p3QuadratureWeights())
#                 Kwn_ss[i,j] = np.asarray([phiI_s(rn,sn)*phiJ_s(rn,sn) for (rn,sn) in rsn]).dot(p3QuadratureWeights())
#                 np.savez("p3_triangle_Krr_Krs_Ksr_Kss.npz",Kwn_rr=Kwn_rr,Kwn_rs=Kwn_rs,Kwn_sr=Kwn_sr,Kwn_ss=Kwn_ss)
            
#     return (Kwn_rr,Kwn_rs,Kwn_sr,Kwn_ss)

def buildVisualizationP3(nodes_x,nodes_y,nodes_z):
    
    # write .node file for TetGen
    node_file = ["# Node count, 3 dim, no attribute, no boundary marker",
                 "%d  3  0  0" % len(nodes_x),
                 "# Node index, node coordinates"]

    for (i,(x,y,z)) in enumerate(zip(nodes_x,nodes_y,nodes_z)):
        node_file.append("%d  %f %f %f" % (i+1,x,y,z))

    facet_list = ["# Part 2 - facet list",
                  "# facet count",
                  "4  0",
                  "# facets",
                  "1",
                  "3  47 48 49",
                  "1",              
                  "3  47 48 50",
                  "1",
                  "3  48 50 49",
                  "1",
                  "3  47 49 50",
                  "",
                  "# Part 3 - hole list",
                  "0",
                  "# Part 4 - no regions",
                  "0"]

    with open("tetP3.poly","w") as f:
        f.write("\n".join(node_file) + "\n" + "\n".join(facet_list))

    with open("tetP3.node","w") as f:
        f.write("\n".join(node_file))



def genP3code():
    autocode = CCodeGen(project="Salvus")
    r,s,t = symbols("r,s,t")
    rstn = p3ReferenceNodes()
    (rn,sn,tn) = rstn
    Np = len(rn)
    phi = p3ElementPolynomials(rstn,r,s,t)
    (dphi_dr,dphi_ds,dphi_dt) = p3ElementDerivatives(phi,r,s,t)
    dphi_dr_rstn = np.zeros((Np,Np))
    dphi_ds_rstn = np.zeros((Np,Np))
    dphi_dt_rstn = np.zeros((Np,Np))
    for i in range(0,Np):
        dphi_dr_i = lambdify((r,s,t),dphi_dr[i])
        dphi_ds_i = lambdify((r,s,t),dphi_ds[i])
        dphi_dt_i = lambdify((r,s,t),dphi_dt[i])        
        for (n,(ri,si,ti)) in enumerate(zip(rn,sn,tn)):
            dphi_dr_rstn[i,n] = dphi_dr_i(ri,si,ti)
            dphi_ds_rstn[i,n] = dphi_ds_i(ri,si,ti)
            dphi_dt_rstn[i,n] = dphi_dt_i(ri,si,ti)

    dphi_dr_rstn = Matrix(dphi_dr_rstn)
    dphi_ds_rstn = Matrix(dphi_ds_rstn)
    dphi_dt_rstn = Matrix(dphi_dt_rstn)
    
    rn = Matrix([ri for (ri,si,ti) in zip(rn,sn,tn)])
    sn = Matrix([si for (ri,si,ti) in zip(rn,sn,tn)])
    tn = Matrix([ti for (ri,si,ti) in zip(rn,sn,tn)])
    wn = Matrix(p3QuadratureWeights())
    rn_routine = autocode.routine("coordinates_p3_tetrahedra_rn",rn,
                                  argument_sequence=None)
    sn_routine = autocode.routine("coordinates_p3_tetrahedra_sn",sn,
                                  argument_sequence=None)
    tn_routine = autocode.routine("coordinates_p3_tetrahedra_tn",tn,
                                  argument_sequence=None)
    wn_routine  = autocode.routine("quadrature_weights_p3_tetrahedra",wn,
                                   argument_sequence=None)

    # NOTE: writes out in Row Major (so we write the
    # transpose). Salvus (via Eigen) assumes Column Major
    dphi_dr_rstn_routine  = autocode.routine("dphi_dr_rstn_p3_tetrahedra",
                                             dphi_dr_rstn.transpose(),
                                             argument_sequence=None)
    dphi_ds_rstn_routine  = autocode.routine("dphi_ds_rstn_p3_tetrahedra",
                                             dphi_ds_rstn.transpose(),
                                             argument_sequence=None)
    dphi_dt_rstn_routine  = autocode.routine("dphi_dt_rstn_p3_tetrahedra",
                                             dphi_dt_rstn.transpose(),
                                             argument_sequence=None)
    
    routines = [rn_routine,
                sn_routine,
                tn_routine,
                wn_routine,
                dphi_dr_rstn_routine,
                dphi_ds_rstn_routine,
                dphi_dt_rstn_routine]

    autocode.write(routines,"p3_tetrahedra",to_files=True)

    buildVisualizationP3(rn,sn,tn)
