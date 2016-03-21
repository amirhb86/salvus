import sympy as sym

# reference quad r=[0,2], s=[0,2] maps into space (x,z) governed by vertices v1,v2,v3,v4

r,s = sym.symbols("r,s")
# r,s->x
v1x,v2x,v3x,v4x = sym.symbols("v1x,v2x,v3x,v4x")
Tx = v1x+(v2x-v1x)*r/2 + (v4x+(v3x-v4x)*r/2 - v1x - (v2x-v1x)*r/2)*s/2
# r,s->z
v1z,v2z,v3z,v4z = sym.symbols("v1z,v2z,v3z,v4z")
Tz = v1z+(v2z-v1z)*r/2 + (v4z+(v3z-v4z)*r/2 - v1z - (v2z-v1z)*r/2)*s/2

dxdr = sym.diff(Tx,r)
dzdr = sym.diff(Tz,r)
dxds = sym.diff(Tx,s)
dzds = sym.diff(Tz,s)

J = sym.Matrix([[dxdr,dzdr],
                [dxds,dzds]])

# Build inverse using 2x2 closed formula
detJxJinv = sym.simplify(sym.Matrix([[dzds,-dzdr],
                                     [-dxds,dxdr]]))

detJ = sym.simplify(J.det())

# detJxJinv = sym.simplify(J.det()*J.inv())
# detJxJinv = Matrix([
# Matrix([
# [ r*(v1z - v2z)/4 + r*(v3z - v4z)/4 - v1z/2 + v4z/2, -s*(v1z - v2z + v3z - v4z)/4 + v1z/2 - v2z/2],
# [-r*(v1x - v2x)/4 - r*(v3x - v4x)/4 + v1x/2 - v4x/2,  s*(v1x - v2x + v3x - v4x)/4 - v1x/2 + v2x/2]])

# detJ = -r*v1x*v3z/8 + r*v1x*v4z/8 + r*v1z*v3x/8 - r*v1z*v4x/8 + r*v2x*v3z/8 - r*v2x*v4z/8 - r*v2z*v3x/8 + r*v2z*v4x/8 - s*v1x*v2z/8 + s*v1x*v3z/8 + s*v1z*v2x/8 - s*v1z*v3x/8 - s*v2x*v4z/8 + s*v2z*v4x/8 + s*v3x*v4z/8 - s*v3z*v4x/8 + v1x*v2z/4 - v1x*v4z/4 - v1z*v2x/4 + v1z*v4x/4 + v2x*v4z/4 - v2z*v4x/4

# Methods for Inverse coordinate transform
Tx = v1x+(v2x-v1x)*(r+1)/2 + (v4x+(v3x-v4x)*(r+1)/2 - v1x - (v2x-v1x)*(r+1)/2)*(s+1)/2
Tz = v1z+(v2z-v1z)*(r+1)/2 + (v4z+(v3z-v4z)*(r+1)/2 - v1z - (v2z-v1z)*(r+1)/2)*(s+1)/2

J = sym.Matrix([[sym.diff(Tx,r), sym.diff(Tx,s)],
                [sym.diff(Tz,r), sym.diff(Tz,s)]])

# J = [-v1x/2 + v2x/2 + (s + 1)*(v1x/2 - v2x/2 + v3x/2 - v4x/2)/2,
# -v1x/2 + v4x/2 - (r + 1)*(-v1x + v2x)/4 + (r + 1)*(v3x - v4x)/4],
# [-v1z/2 + v2z/2 + (s + 1)*(v1z/2 - v2z/2 + v3z/2 - v4z/2)/2, -v1z/2
# + v4z/2 - (r + 1)*(-v1z + v2z)/4 + (r + 1)*(v3z - v4z)/4]])
