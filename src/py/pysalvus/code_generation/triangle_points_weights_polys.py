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

def p3QuadratureWeights():
    ws = 2*( 919*sqrt(7) + 2471 )/(124080*sqrt(7) + 330960)
    wa = 2*( sqrt(7)*(2+sqrt(7))**4 )/( 25280 + 9520*sqrt(7) )
    wb = 2*( 147 + 42*sqrt(7) ) / ( 400*sqrt(7) + 1280 )

    # weights correspond to nodal layout seen below
    quadrature_weights = np.asarray([wb,wb,wb,
                                     wa,wa,wa,wa,wa,wa,
                                     ws,ws,ws])

    return quadrature_weights

def p3ReferenceNodes():
    
    # P3 reference triangle in PETSc ordering (face,edges,vertices)
    #    11
    #    | \--
    #    |    \--
    #    7       6--
    #    |          \--
    #    |             \-
    #    |     2         \--
    #    8                  5--
    #    |     0    1          \--
    #    |                        \--
    #    9------3---------4--------- 10

    # petsc ordering
    rsn = [(-0.5853,-0.5853),(0.1706,-0.5853),(-0.5853,0.1706), # G_{1,2,3}
           (-0.413,-1.0),(0.413,-1.0),(0.413,-0.413),(-0.413,0.413),(-1.0,0.413),(-1.0,-0.413), # M_{(123,123)}
           (-1.0,-1.0),(1.0,-1.0),(-1.0,1.0)] # S_{1,2,3}
    rsn = np.asarray(rsn)
    return rsn

def p3ElementPolynomials(rsn,r,s):
    b = (-r-s)*(s+1)*(r+1)
    P1b = np.asarray([r,s])*b
    P3 = np.asarray([1+0*r,r,s,r*s,r**2,s**2,r**3,s**3,
                     r**2*s,r*s**2])

    P3plus_P1b = np.concatenate((P3,P1b))
    Np=len(P3plus_P1b)
    V = np.zeros((Np,Np))
    for i in range(0,Np):
        (rn,sn) = rsn[i]
        V[i,:] =  np.asarray([P3plus_P1b_j.subs({r:rn,s:sn}) for P3plus_P1b_j in P3plus_P1b])

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
        phi_i = np.dot(P3plus_P1b, p_ai)[0]
        phi.append(phi_i)

    return phi

def p3ElementDerivatives(phi,r,s):
    deriv_r = []
    deriv_s = []
    for phi_i in phi:
        dphi_dr = diff(phi_i,r)
        dphi_ds = diff(phi_i,s)
        deriv_r.append(dphi_dr)
        deriv_s.append(dphi_ds)
    return (deriv_r,deriv_s)

def genK(phi,wn,r,s,rsn):
    Np = len(wn)
    Kwn_rr = np.zeros((Np,Np))
    Kwn_rs = np.zeros((Np,Np))
    Kwn_sr = np.zeros((Np,Np))
    Kwn_ss = np.zeros((Np,Np))
    if os.path.isfile("p3_triangle_Krr_Krs_Ksr_Kss.npz"):
        print("Using precomputed Krr,etc")
        Kall = np.load("Krr_Krs_Ksr_Kss.npz")
        Kwn_rr = Kall["Kwn_rr"]
        Kwn_rs = Kall["Kwn_rs"]
        Kwn_sr = Kall["Kwn_sr"]
        Kwn_ss = Kall["Kwn_ss"]
    else:
        print("Computing Krr,etc and storing them")
        for i in range(0,Np):
            print("i={}".format(i))
            for j in range(0,Np):
                # exact version
                # phiI_x_phiJ_rr = diff(phi[i],r)*diff(phi[j],r)
                # phiI_x_phiJ_ss = diff(phi[i],s)*diff(phi[j],s)
                # Kexact[i,j] = integrate(phiI_x_phiJ_rr,(r,-1,-s),(s,-1,1)) + integrate(phiI_x_phiJ_ss,(r,-1,-s),(s,-1,1))
                # phiI_x_phiJ_rr = diff(phi[i],r)*diff(phi[j],r)
                # phiI_x_phiJ_ss = diff(phi[i],s)*diff(phi[j],s)
                # phiI_x_phiJ_rs = diff(phi[i],r)*diff(phi[j],s)
                # phiI_x_phiJ_sr = diff(phi[i],s)*diff(phi[j],r)
                # Kwn_rr[i,j] = integrate(phiI_x_phiJ_rr,(r,-1,-s),(s,-1,1))
                # Kwn_rs[i,j] = integrate(phiI_x_phiJ_rs,(r,-1,-s),(s,-1,1))
                # Kwn_sr[i,j] = integrate(phiI_x_phiJ_sr,(r,-1,-s),(s,-1,1))
                # Kwn_ss[i,j] = integrate(phiI_x_phiJ_ss,(r,-1,-s),(s,-1,1))

                # version using quadrature
                phiI_r = lambdify((r,s),diff(phi[i],r))
                phiI_s = lambdify((r,s),diff(phi[i],s))
                phiJ_r = lambdify((r,s),diff(phi[j],r))
                phiJ_s = lambdify((r,s),diff(phi[j],s))
                Kwn_rr[i,j] = np.asarray([phiI_r(rn,sn)*phiJ_r(rn,sn) for (rn,sn) in rsn]).dot(p3QuadratureWeights())
                Kwn_rs[i,j] = np.asarray([phiI_r(rn,sn)*phiJ_s(rn,sn) for (rn,sn) in rsn]).dot(p3QuadratureWeights())
                Kwn_sr[i,j] = np.asarray([phiI_s(rn,sn)*phiJ_r(rn,sn) for (rn,sn) in rsn]).dot(p3QuadratureWeights())
                Kwn_ss[i,j] = np.asarray([phiI_s(rn,sn)*phiJ_s(rn,sn) for (rn,sn) in rsn]).dot(p3QuadratureWeights())
                np.savez("p3_triangle_Krr_Krs_Ksr_Kss.npz",Kwn_rr=Kwn_rr,Kwn_rs=Kwn_rs,Kwn_sr=Kwn_sr,Kwn_ss=Kwn_ss)
            
    return (Kwn_rr,Kwn_rs,Kwn_sr,Kwn_ss)


def genP3code():
    autocode = CCodeGen(project="Salvus")
    r,s = symbols("r,s")
    rsn = p3ReferenceNodes()
    Np = len(rsn)
    phi = p3ElementPolynomials(rsn,r,s)
    (dphi_dr,dphi_ds) = p3ElementDerivatives(phi,r,s)
    dphi_dr_rsn = np.zeros((Np,Np))
    dphi_ds_rsn = np.zeros((Np,Np))
    for i in range(0,Np):
        dphi_dr_i = lambdify((r,s),dphi_dr[i])
        dphi_ds_i = lambdify((r,s),dphi_ds[i])
        for (n,(rn,sn)) in enumerate(rsn):
            dphi_dr_rsn[i,n] = dphi_dr_i(rn,sn)
            dphi_ds_rsn[i,n] = dphi_ds_i(rn,sn)

    dphi_dr_rsn = Matrix(dphi_dr_rsn)
    dphi_ds_rsn = Matrix(dphi_ds_rsn)
    print(dphi_dr_rsn[0,:])
    rn = Matrix([ri for (ri,si) in rsn])
    sn = Matrix([si for (ri,si) in rsn])
    wn = Matrix(p3QuadratureWeights())
    rn_routine = autocode.routine("coordinates_p3_triangle_rn",rn,
                                  argument_sequence=None)
    sn_routine = autocode.routine("coordinates_p3_triangle_sn",sn,
                                  argument_sequence=None)
    wn_routine  = autocode.routine("quadrature_weights_p3_triangle",wn,
                                   argument_sequence=None)

    # NOTE: by default, it writes out in Row Major (so we write the
    # transpose). Salvus (via Eigen) assumes Column Major
    dphi_dr_rsn_routine  = autocode.routine("dphi_dr_rsn_p3_triangle",
                                            dphi_dr_rsn.transpose(),
                                            argument_sequence=None)
    dphi_ds_rsn_routine  = autocode.routine("dphi_ds_rsn_p3_triangle",
                                            dphi_ds_rsn.transpose(),
                                            argument_sequence=None)

    routines = [rn_routine,
                sn_routine,
                wn_routine,
                dphi_dr_rsn_routine,
                dphi_ds_rsn_routine]

    autocode.write(routines,"p3_triangle",to_files=True)
