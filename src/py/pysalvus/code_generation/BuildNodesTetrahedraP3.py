import numpy as np
import pylab as pl
import numpy.linalg

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

def BuildNodesTetrahedraP3(plot):
    
    def plotTet(ax):

        x_1 = [0,1,0,0]
        y_1 = [0,0,1,0]
        z_1 = [0,0,0,0]

        x_2 = [0,1,0,0]
        y_2 = [0,0,0,0]
        z_2 = [0,0,1,0]

        x_3 = [0,0,0,0]
        y_3 = [0,1,0,0]
        z_3 = [0,0,1,0]

        x_4 = [1,0,0,1]
        y_4 = [0,1,0,0]
        z_4 = [0,0,1,0]

        tet_x = x_1+x_2+x_3+x_4
        tet_y = y_1+y_2+y_3+y_4
        tet_z = z_1+z_2+z_3+z_4

        ax.plot(tet_x,tet_y,tet_z)

        # midpoint lines
        # side 1 x-y
        ax.plot([0,1/2.],[0,1/2.],linestyle="dashed",color="black")
        ax.plot([1,0],[0,1/2.],linestyle="dashed",color="black")
        ax.plot([0,1/2.],[1,0],linestyle="dashed",color="black")

        # side 2 x-z
        ax.plot([0,1/2.],[0,0],[0,1/2.],linestyle="dashed",color="black")
        ax.plot([1,0],[0,0],[0,1/2.],linestyle="dashed",color="black")
        ax.plot([0,1/2.],[0,0],[1,0],linestyle="dashed",color="black")

        # side 3 z-y
        ax.plot([0,0],[0,1/2.],[0,1/2.],linestyle="dashed",color="black")
        ax.plot([0,0],[0,1/2.],[1,0],linestyle="dashed",color="black")
        ax.plot([0,0],[1,0],[0,1/2.],linestyle="dashed",color="black")

        # side 4 x-y-z
        ax.plot([1,0],[0,1/2.],[0,1/2.],linestyle="dashed",color="black")
        ax.plot([0,1/2.],[0,1/2.],[1,0],linestyle="dashed",color="black")
        ax.plot([0,1/2.],[1,0],[0,1/2.],linestyle="dashed",color="black")

    def makeT((p1,p2,p3)):
        return np.asmatrix([[p1[0]-p3[0],p2[0]-p3[0]],[p1[1]-p3[1],p2[1]-p3[1]]])

    def toBC(invT,r,r3):
        l12 = invT.dot(r-r3)    
        return (l12[0,0],l12[0,1],1-l12[0,0]-l12[0,1])

    def fromBC(T,r3,l1,l2):
        r_m_r3 = T.dot(np.asarray([l1,l2]))
        return r_m_r3 + r3


    def faceCoords_bc(beta,r1,r2,r3):
        (x1,y1) = r1
        (x2,y2) = r2
        (x3,y3) = r3
        T = makeT((r1,r2,r3))    
        # p1
        (l1,l2,l3) = toBC(np.linalg.inv(T),np.asarray([beta,beta]),np.asarray([x3,y3]))
        p1bc = (l1,l2,l3)
        beta_bc = l1
        p2bc = ((1-beta_bc)/2,beta_bc,(1-beta_bc)/2)
        p3bc = ((1-beta_bc)/2,(1-beta_bc)/2,beta_bc)
        p1 = fromBC(T,np.asarray([x3,y3]),p1bc[0],p1bc[1])
        p2 = fromBC(T,np.asarray([x3,y3]),p2bc[0],p2bc[1])
        p3 = fromBC(T,np.asarray([x3,y3]),p3bc[0],p3bc[1])
        return (np.asarray(p1)[0],np.asarray(p2)[0],np.asarray(p3)[0])
    
    def edgeCoords(alpha,r1,r2):

        p1 = (1-alpha)*r1 + alpha*r2
        p2 = alpha*r1 + (1-alpha)*r2

        return (p1,p2)

    def faceCoords(beta,r1,r2,r3):

        r1m = (r2+r3)/2.0
        r2m = (r1+r3)/2.0
        r3m = (r1+r2)/2.0
        p1 = (1-2*beta)*r1 + (2*beta)*r1m
        p2 = (1-2*beta)*r2 + (2*beta)*r2m
        p3 = (1-2*beta)*r3 + (2*beta)*r3m

        return (p1,p2,p3)

    def getFace(r1,r2,r3,beta,ax,plot=True):

        (p1a,p1b,p1c) = faceCoords(beta,r1,r2,r3)

        return ([p1a[0],p1b[0],p1c[0]],[p1a[1],p1b[1],p1c[1]],[p1a[2],p1b[2],p1c[2]])

    def getEdge(r1,r2,alpha,ax,plot=True):

        (p1a,p1b) = edgeCoords(alpha,r1,r2)
        
        return ([p1a[0],p1b[0]],[p1a[1],p1b[1]],[p1a[2],p1b[2]])

    # interior points
    def intCoordsA(gamma,r1,r2,r3,r4,ax,plot=True):

        r1m = (r2+r3+r4)/3.0
        r2m = (r1+r3+r4)/3.0
        r3m = (r1+r2+r4)/3.0
        r4m = (r1+r2+r3)/3.0

        p1 = (1-3*gamma)*r1 + (3*gamma)*r1m
        p2 = (1-3*gamma)*r2 + (3*gamma)*r2m
        p3 = (1-3*gamma)*r3 + (3*gamma)*r3m
        p4 = (1-3*gamma)*r4 + (3*gamma)*r4m
        
        return (p1,p2,p3,p4)

    def getInteriorA(gamma,r1,r2,r3,r4,ax,plot=True):

        (pA,pB,pC,pD) = intCoordsA(gamma,r1,r2,r3,r4,ax,plot)
        
        return ([pA[0],pB[0],pC[0],pD[0]],[pA[1],pB[1],pC[1],pD[1]],[pA[2],pB[2],pC[2],pD[2]])

    def intCoordsB(delta,r1,r2,r3,r4,ax,plot=True):

        r14m = (r1+r4)/2.0
        r12m = (r1+r2)/2.0
        r13m = (r1+r3)/2.0
        r24m = (r2+r4)/2.0
        r23m = (r2+r3)/2.0
        r34m = (r3+r4)/2.0

        p1 = (1-2*delta)*r14m + (2*delta)*r23m
        p2 = (1-2*delta)*r12m + (2*delta)*r34m
        p3 = (1-2*delta)*r13m + (2*delta)*r24m

        p4 = (1-2*delta)*r23m + (2*delta)*r14m
        p5 = (1-2*delta)*r34m + (2*delta)*r12m
        p6 = (1-2*delta)*r24m + (2*delta)*r13m
        
        return (p1,p2,p3,p4,p5,p6)


    def getInteriorB(delta,r1,r2,r3,r4,ax,plot=True):

        (p1,p2,p3,p4,p5,p6) = intCoordsB(delta,r1,r2,r3,r4,ax,plot)
        # if plot:
        #     ax.plot([p1[0]],[p1[1]],[p1[2]],marker="o",color="yellow")
        #     ax.plot([p2[0]],[p2[1]],[p2[2]],marker="o",color="yellow")
        #     ax.plot([p3[0]],[p3[1]],[p3[2]],marker="o",color="yellow")
        #     ax.plot([p4[0]],[p4[1]],[p4[2]],marker="o",color="yellow")
        #     ax.plot([p5[0]],[p5[1]],[p5[2]],marker="o",color="yellow")
        #     ax.plot([p6[0]],[p6[1]],[p6[2]],marker="o",color="yellow")
        
        return ([p1[0],p2[0],p3[0],p4[0],p5[0],p6[0]],[p1[1],p2[1],p3[1],p4[1],p5[1],p6[1]],[p1[2],p2[2],p3[2],p4[2],p5[2],p6[2]])



    fig = plt.figure(1)
    ax = fig.gca(projection='3d')

    beta1 = 0.1972862280257976
    beta2 = 0.4256461243139345

    # side 1 x-y
    r1 = np.asarray([0,0,0])
    r2 = np.asarray([0,1,0])
    r3 = np.asarray([1,0,0])
    (x1a,y1a,z1a) = getFace(r1,r2,r3,beta1,ax,plot)
    (x1b,y1b,z1b) = getFace(r1,r2,r3,beta2,ax,plot)
    
    # side 2 z-y
    r1 = np.asarray([0,0,0])
    r2 = np.asarray([0,0,1])
    r3 = np.asarray([0,1,0])
    (x2a,y2a,z2a) = getFace(r1,r2,r3,beta1,ax,plot)
    (x2b,y2b,z2b) = getFace(r1,r2,r3,beta2,ax,plot)

    # side 3 x-z
    r1 = np.asarray([0,0,0])
    r2 = np.asarray([1,0,0])
    r3 = np.asarray([0,0,1])
    (x3a,y3a,z3a) = getFace(r1,r2,r3,beta1,ax,plot)
    (x3b,y3b,z3b) = getFace(r1,r2,r3,beta2,ax,plot)
    
    # side 4 x-y-z
    r1 = np.asarray([1,0,0])
    r2 = np.asarray([0,1,0])
    r3 = np.asarray([0,0,1])
    (x4a,y4a,z4a) = getFace(r1,r2,r3,beta1,ax,plot)
    (x4b,y4b,z4b) = getFace(r1,r2,r3,beta2,ax,plot)

    # face construction complete, assemble face nodes
    facepoints_x = x1a+x1b+x2a+x2b+x3a+x3b+x4a+x4b
    facepoints_y = y1a+y1b+y2a+y2b+y3a+y3b+y4a+y4b
    facepoints_z = z1a+z1b+z2a+z2b+z3a+z3b+z4a+z4b

    # edges
    alpha = 0.2928294047674109

    r1 = np.asarray([0,0,0])
    r2 = np.asarray([0,1,0])
    (x1,y1,z1)= getEdge(r1,r2,alpha,ax,plot)

    r1 = np.asarray([0,1,0])
    r2 = np.asarray([1,0,0])
    (x2,y2,z2)= getEdge(r1,r2,alpha,ax,plot)

    r1 = np.asarray([1,0,0])
    r2 = np.asarray([0,0,0])
    (x3,y3,z3)= getEdge(r1,r2,alpha,ax,plot)

    r1 = np.asarray([0,0,0])
    r2 = np.asarray([0,0,1])
    (x4,y4,z4)= getEdge(r1,r2,alpha,ax,plot)
    
    r1 = np.asarray([0,0,1])
    r2 = np.asarray([0,1,0])
    (x5,y5,z5)= getEdge(r1,r2,alpha,ax,plot)
    
    r1 = np.asarray([1,0,0])
    r2 = np.asarray([0,0,1])
    (x6,y6,z6)= getEdge(r1,r2,alpha,ax,plot)
    
    edgepoints_x = x1+x2+x3+x4+x5+x6
    edgepoints_y = y1+y2+y3+y4+y5+y6
    edgepoints_z = z1+z2+z3+z4+z5+z6

    corners_x = [0,0,1,0]
    corners_y = [0,1,0,0]
    corners_z = [0,0,0,1]

    # corner to face midline
    r1 = np.asarray([0,0,0])
    r2 = np.asarray([1,0,0])
    r3 = np.asarray([0,1,0])
    r4 = np.asarray([0,0,1])

    gamma = 0.9503775858394107e-01
    (intA_x,intA_y,intA_z) = getInteriorA(gamma,r1,r2,r3,r4,ax,plot)

    delta = 0.1252462362578136
    (intB_x,intB_y,intB_z) = getInteriorB(delta,r1,r2,r3,r4,ax,plot)

    int_x = intA_x + intB_x
    int_y = intA_y + intB_y
    int_z = intA_z + intB_z

    # pts given in order:
    # {cell, faces, edges, vertices}
    # {interior, faces, edges, vertices}
    
    all_points_x = int_x + facepoints_x + edgepoints_x + corners_x
    all_points_y = int_y + facepoints_y + edgepoints_y + corners_y
    all_points_z = int_z + facepoints_z + edgepoints_z + corners_z

    if plot:
        
        # plot tetrahedral frame
        plotTet(ax)
        for (i,(x,y,z)) in enumerate(zip(all_points_x,all_points_y,all_points_z)):
            if i>=34 and i < 50:
                ax.text(x,y,z,"{}".format(i))
            # ax.plot(all_points_x,all_points_y,all_points_z,marker="o",color="black",linestyle='None')
        plt.xlabel("x")
        plt.ylabel("y")
        # plt.zlabel("z")
        plt.show()

    np.savez("RST_p3tetrahedra.npz",r=all_points_x,s=all_points_y,t=all_points_z)
        
    return (all_points_x,all_points_y,all_points_z)
