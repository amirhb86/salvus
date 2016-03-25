import sympy as sym

# inverse coordinate transform

# from wikipedia: https://en.wikipedia.org/wiki/Barycentric_coordinate_system
# Given triangle with vertices (x1,z1)...

x1, x2, x3, z1, z2, z3, x, z = sym.symbols("x1 x2 x3 z1 z2 z3 x z")

# for reference triangle
Tref = sym.Matrix([[0, 2],
                [-2,-2]])

# Tinv = T.inv()

# for any triangle (x1,y1),...
T = sym.Matrix([[x1-x3,x2-x3],
                [z1-z3,z2-z3]])

# barycentric coordinates l1,l2,l3 (l1+l2+l3=1) represent any point r_n = l1*v1+l2*v2+l3*v3 given any vertices.
# T*[l1;l2] = xy - v3, where xy is the desired point, and v3 is the vertex 3
# thus l12 = T**(-1)*xy - T**(-1)*v3

# Solving for l1,l2,l3 will allow us to get reference coords (r,s) = l1*v1_ref*l2*v2_ref*l3*v3_ref

l12 = T.inv()*sym.Matrix([x,z]) - T.inv()*sym.Matrix([x3,z3])
l1 = l12[0]
l2 = l12[1]
l3 = 1-l1-l2

# multiply barycentric coords l1,l2,l3 with reference vertices to get corresponding reference point.
rn = sym.simplify(l1*(-1) + l2*(1) * l3*(-1))
sn = sym.simplify(l1*(-1) + l2*(-1) * l3*(1))

print("r={},\ns={}".format(rn,sn))

#r=((-x + x3)*(z2 - z3)*((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)) + (x2 - x3)*(z - z3)*((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)) + (x*(z1 - z3) - x3*(z1 - z3) - z*(x1 - x3) + z3*(x1 - x3))*(x*(z1 - z3) - x3*(z1 - z3) - z*(x1 - x3) + z*(x2 - x3) + z3*(x1 - x3) - z3*(x2 - x3) + (-x + x3)*(z2 - z3) + (x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)))/((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3))**2,
#s=((-x + x3)*(z2 - z3)*((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)) + (x2 - x3)*(z - z3)*((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)) + (x*(z1 - z3) - x3*(z1 - z3) - z*(x1 - x3) + z3*(x1 - x3))*(x*(z1 - z3) - x3*(z1 - z3) - z*(x1 - x3) + z*(x2 - x3) + z3*(x1 - x3) - z3*(x2 - x3) + (-x + x3)*(z2 - z3) + (x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)))/((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3))**2

# transform from reference point (r,s) to triangle point (x,z)
r,s = sym.symbols("r,s")
# l12 = Tinv*(r;s)-Tinv*(v3)
l12 = Tref.inv()*sym.Matrix([r,s])-Tref.inv()*sym.Matrix([-1,1]) # (-1,1) is v3 in reference triangle

l1 = sym.simplify(l12[0])
l2 = sym.simplify(l12[1])
l3 = sym.simplify(1-l12[0]-l12[1])

xn = sym.simplify(l1*x1 + l2*x2 + l3*x3)
zn = sym.simplify(l1*z1 + l2*z2 + l3*z3)

print("xn={},\nzn={}".format(xn,zn))
# l1=-r/2 - s/2,
# l2=r/2 + 1/2,
# l3=s/2 + 1/2
# xn=-x1*(r + s)/2 + x2*(r + 1)/2 + x3*(s + 1)/2,
# zn=-z1*(r + s)/2 + z2*(r + 1)/2 + z3*(s + 1)/2
