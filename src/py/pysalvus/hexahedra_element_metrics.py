import sympy as sym

# reference hex, with r,s,t=[-1,1]x[-1,1]x[-1,1]
#
#                                  (v7)                    (v6)
#                                    /--------------d------/+
#                                 /--|                  /-- |
#                             /---   |              /---    |
#                          /--       |       f   /--        |
#                       /--          |         --           |
#                (v4) +---------------b--------+ (v5)       |
#                     |              |         |            |
#                     |              |         |            |
#                     |              |       g |            |
#   ^                 |              |         |            |
#   | (t)             |              |         |            |
#   |                 |            --+-------- |----c-------+ (v2)
#   |                 |        ---/ (v1)       |         /--
#   |                 |     --/              e |     /---
#   |         (s)     | ---/                   |  /--
#   |         /-      +/--------------a--------+--
#   |     /---       (v0)                    (v3)
#   | /---
#   +-----------------> (r)
#
#

r,s,t = sym.symbols("r,s,t")
v0,v1,v2,v3,v4,v5,v6,v7 = sym.symbols("v0,v1,v2,v3,v4,v5,v6,v7")


# note (r+1)/2.0 to bring range from [-1,1] to [0,1]
a = v0 + (v3-v0)*(r+1.0)/2.0
b = v4 + (v5-v4)*(r+1.0)/2.0
c = v1 + (v2-v1)*(r+1.0)/2.0
d = v7 + (v6-v7)*(r+1.0)/2.0

e = a + (c-a)*(s+1.0)/2.0
f = b + (d-b)*(s+1.0)/2.0

g = e + (f-e)*(t+1.0)/2.0

print("mapping g={}".format(sym.simplify(g)))

vx = [0,0,1,1,1.0/2.0,3.0/2.0,3.0/2.0,1.0/2.0]
vy = [0,1,1,0,1.0/2.0,1.0/2.0,3.0/2.0,3.0/2.0]
vz = [0,0,0,0,1.0/2.0,1.0/2.0,1.0/2.0,1.0/2.0]

# T(-1.000003,1.000000,0.447214) -> (0.276393,-0.000000,-0.000000) vs. goal (0.276393,0.000000,0.000000)
# (-0.447214,-1.000000,-1.000000) =?= (-1.000003,1.000000,0.447214)

# recent
# T(-0.447214,1.000000,-1.000000) -> (0.723607,0.000000,0.276393)

# g=v0 - 0.5*(r + 1.0)*(v0 - v3) - 0.5*(s + 1.0)*(v0 - v1 - 0.5*(r + 1.0)*(v0 - v3) + 0.5*(r + 1.0)*(v1 - v2)) + 0.5*(t + 1.0)*(-v0 + v4 + 0.5*(r + 1.0)*(v0 - v3) - 0.5*(r + 1.0)*(v4 - v5) + 0.5*(s + 1.0)*(v0 - v1 - 0.5*(r + 1.0)*(v0 - v3) + 0.5*(r + 1.0)*(v1 - v2)) + 0.5*(s + 1.0)*(-v4 + v7 + 0.5*(r + 1.0)*(v4 - v5) + 0.5*(r + 1.0)*(v6 - v7)))

gx = g.subs({v0:vx[0],v1:vx[1],v2:vx[2],v3:vx[3],v4:vx[4],v5:vx[5],v6:vx[6],v7:vx[7]})
gy = g.subs({v0:vy[0],v1:vy[1],v2:vy[2],v3:vy[3],v4:vy[4],v5:vy[5],v6:vy[6],v7:vy[7]})
gz = g.subs({v0:vz[0],v1:vz[1],v2:vz[2],v3:vz[3],v4:vz[4],v5:vz[5],v6:vz[6],v7:vz[7]})



## Derive jacobian terms

# `g` maps us from reference coordinate system into physical coordinates based on input vertices `v0-v7`.

# use x,y,z vertex names to make copy-paste easier.
v0x,v1x,v2x,v3x,v4x,v5x,v6x,v7x = sym.symbols("v0x,v1x,v2x,v3x,v4x,v5x,v6x,v7x")
Tx = g.subs({v0:v0x,v1:v1x,v2:v2x,v3:v3x,v4:v4x,v5:v5x,v6:v6x,v7:v7x})
v0y,v1y,v2y,v3y,v4y,v5y,v6y,v7y = sym.symbols("v0y,v1y,v2y,v3y,v4y,v5y,v6y,v7y")
Ty = g.subs({v0:v0y,v1:v1y,v2:v2y,v3:v3y,v4:v4y,v5:v5y,v6:v6y,v7:v7y})
v0z,v1z,v2z,v3z,v4z,v5z,v6z,v7z = sym.symbols("v0z,v1z,v2z,v3z,v4z,v5z,v6z,v7z")
Tz = g.subs({v0:v0z,v1:v1z,v2:v2z,v3:v3z,v4:v4z,v5:v5z,v6:v6z,v7:v7z})

dxdr = sym.diff(Tx,r)
dydr = sym.diff(Ty,r)
dzdr = sym.diff(Tz,r)
dxds = sym.diff(Tx,s)
dyds = sym.diff(Ty,s)
dzds = sym.diff(Tz,s)
dxdt = sym.diff(Tx,t)
dydt = sym.diff(Ty,t)
dzdt = sym.diff(Tz,t)

# J = sym.simplify(sym.Matrix([[dxdr, dydr, dzdr],
#                              [dxds, dyds, dzds],
#                              [dxdt, dydt, dzdt]]))

# print("J={}".format(J))

# # copy-pasta for C++
# print("PetscReal dxdr={};".format(dxdr))
# print("PetscReal dydr={};".format(dydr))
# print("PetscReal dzdr={};".format(dzdr))
# print("PetscReal dxds={};".format(dxds))
# print("PetscReal dyds={};".format(dyds))
# print("PetscReal dzds={};".format(dzds))
# print("PetscReal dxdt={};".format(dxdt))
# print("PetscReal dydt={};".format(dydt))
# print("PetscReal dzdt={};".format(dzdt))

# detJ = J.det()
# # print("detJ={}".format(detJ))
# invJ = J.inv()
# # print("invJ={}".format(detJ))

# detJxJinv = sym.simplify(detJ*invJ)
# # print("detJxInvJ={}".format(detJxInvJ))


# inverse Jacobian
dxdr_ = sym.symbols("dxdr")
dydr_ = sym.symbols("dydr")
dzdr_ = sym.symbols("dzdr")
dxds_ = sym.symbols("dxds")
dyds_ = sym.symbols("dyds")
dzds_ = sym.symbols("dzds")
dxdt_ = sym.symbols("dxdt")
dydt_ = sym.symbols("dydt")
dzdt_ = sym.symbols("dzdt")

J_ = sym.Matrix([
    [dxdr_, dydr_, dzdr_],
    [dxds_, dyds_, dzds_],
    [dxdt_, dydt_, dzdt_]])

# detJxInvJ_ = sym.simplify(J_.det()*J_.inv())


# linear interpolation for velocity
# reference hex, with r,s,t=[-1,1]x[-1,1]x[-1,1]
#
#                                               (v7)                    (v6)
#                                                 /--------------d------/+
#                                              /--|                  /-- |
#                                          /---   |              /---    |
#                                       /--       |       f   /--        |
#                                    /--          |         --           |
#                             (v4) +---------------b--------+ (v5)       |        \gamma
#                                  |              |         |            |          ^
#                                  |              |         |            |          |
#                                  |              |       g |            |          |
#   ^                              |              |         |            |          |
#   | (t)                          |              |         |            |          |
#   |                              |            --+-------- |----c-------+ (v2)     |
#   |                              |        ---/ (v1)       |         /--          ---
#   |                              |     --/              e |     /---               /> \beta
#   |         (s)                  | ---/                   |  /--                 /-
#   |         /-                   +/--------------a--------+--                  --
#   |     /---                    (v0)                    (v3)
#   | /---                         |--------------> \alpha
#   +-----------------> (r)
#
#

r,s,t,v0,v1,v2,v3,v4,v5,v6,v7 = sym.symbols("r,s,t,v0,v1,v2,v3,v4,v5,v6,v7")
# alpha,beta,gamma are in range [0,1]
alpha = (r+1)/2.0
beta = (s+1)/2.0
gamma = (t+1)/2.0

a = alpha*v3 + (1-alpha)*v0
b = alpha*v5 + (1-alpha)*v4
c = alpha*v2 + (1-alpha)*v1
d = alpha*v6 + (1-alpha)*v7

e = beta*c + (1-beta)*a
f = beta*d + (1-beta)*b

g = gamma*f + (1-gamma)*e

g_full = sym.collect(sym.expand(g),[v0,v1,v2,v3,v4,v5,v6,v7])
print("g={}".format(g_full))
g_terms = sym.collect(sym.expand(g),[v0,v1,v2,v3,v4,v5,v6,v7],evaluate=False)

v0t = g_terms[v0]
v1t = g_terms[v1]
v2t = g_terms[v2]
v3t = g_terms[v3]
v4t = g_terms[v4]
v5t = g_terms[v5]
v6t = g_terms[v6]
v7t = g_terms[v7]
print("vector=[{},{},{},{},{},{},{},{}]".format(v0t,v1t,v2t,v3t,v4t,v5t,v6t,v7t))
# vector=[-0.125*r*s*t + 0.125*r*s + 0.125*r*t - 0.125*r + 0.125*s*t - 0.125*s - 0.125*t + 0.125,0.125*r*s*t - 0.125*r*s + 0.125*r*t - 0.125*r - 0.125*s*t + 0.125*s - 0.125*t + 0.125,-0.125*r*s*t + 0.125*r*s - 0.125*r*t + 0.125*r - 0.125*s*t + 0.125*s - 0.125*t + 0.125,0.125*r*s*t - 0.125*r*s - 0.125*r*t + 0.125*r + 0.125*s*t - 0.125*s - 0.125*t + 0.125,0.125*r*s*t + 0.125*r*s - 0.125*r*t - 0.125*r - 0.125*s*t - 0.125*s + 0.125*t + 0.125,-0.125*r*s*t - 0.125*r*s + 0.125*r*t + 0.125*r - 0.125*s*t - 0.125*s + 0.125*t + 0.125,0.125*r*s*t + 0.125*r*s + 0.125*r*t + 0.125*r + 0.125*s*t + 0.125*s + 0.125*t + 0.125,-0.125*r*s*t - 0.125*r*s - 0.125*r*t - 0.125*r + 0.125*s*t + 0.125*s + 0.125*t + 0.125]



