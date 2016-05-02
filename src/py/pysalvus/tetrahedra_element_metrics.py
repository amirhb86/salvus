import sympy as sym

# inverse coordinate transform

# from wikipedia: https://en.wikipedia.org/wiki/Barycentric_coordinate_system
# Given triangle with vertices (x1,z1)...

x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, x, y, z = sym.symbols("x1 x2 x3 x4 y1 y2 y3 y4 z1 z2 z3 z4 x y z")

# for any triangle (x1,y1,z1),...
T = sym.Matrix([[x1-x4,x2-x4,x3-x4],
                [y1-y4,y2-y4,y3-y4],
                [z1-z4,z2-z4,z3-z4]])
                
Tref = T.subs({x1:0,y1:0,z1:0,x2:0,y2:1,z2:0,x3:1,y3:0,z3:0,x4:0,y4:0,z4:1})

# barycentric coordinates l1,l2,l3 (l1+l2+l3=1) represent any point r_n = l1*v1+l2*v2+l3*v3 given any vertices.
# T*[l1;l2] = xy - v3, where xy is the desired point, and v3 is the vertex 3
# thus l12 = T**(-1)*xy - T**(-1)*v3

# Solving for l1,l2,l3 will allow us to get reference coords (r,s) = l1*v1_ref*l2*v2_ref*l3*v3_ref

l123 = T.inv()*sym.Matrix([x,y,z]) - T.inv()*sym.Matrix([x4,y4,z4])
l1 = sym.simplify(l123[0])
l2 = sym.simplify(l123[1])
l3 = sym.simplify(l123[2])
l4 = 1-l1-l2-l3

print("l1={},\nl2={},\nl3={},\nl4={}".format(l1,l2,l3,l4))

# multiply barycentric coords l1,l2,l3 with reference vertices to get corresponding reference point.
rn = sym.simplify(l3) # l1*(0) + l2*(0) * l3*(1) + l4*(0)
sn = sym.simplify(l2) # l1*(0) + l2*(1) * l3*(0) + l4*(0)
tn = sym.simplify(l4) # l1*(0) + l2*(0) * l3*(0) + l4*(1)

print("auto r={};\nauto s={};\nauto t={};".format(rn,sn,tn))

# transform from reference point (r,s) to triangle point (x,z)
r,s,t = sym.symbols("r,s,t")
# l123 = Tinv*(r;s;t)-Tinv*(v4)
l123 = Tref.inv()*sym.Matrix([r,s,t])-Tref.inv()*sym.Matrix([0,0,1]) # (0,0,1) is v4 in reference tet

l1 = sym.simplify(l123[0])
l2 = sym.simplify(l123[1])
l3 = sym.simplify(l123[2])
l4 = sym.simplify(1-l1-l2-l3)

xn = sym.simplify(l1*x1 + l2*x2 + l3*x3 + l4*x4)
yn = sym.simplify(l1*y1 + l2*y2 + l3*y3 + l4*y4)
zn = sym.simplify(l1*z1 + l2*z2 + l3*z3 + l4*z4)

print("xn={},\nyn={},\nzn={}".format(xn,yn,zn))

# vertex interpolation
v0,v1,v2,v3 = sym.symbols("v0,v1,v2,v3")
vector = [l1,l2,l3,l4]
print("interp_vector={}".format(vector))
