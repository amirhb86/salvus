
import sympy as sym

eps, eta, x, z = sym.symbols('epsilon eta x z')
n0 = (1.0/4.0) * (1.0 - eps) * (1.0 - eta)
n1 = (1.0/4.0) * (1.0 + eps) * (1.0 - eta)
n2 = (1.0/4.0) * (1.0 - eps) * (1.0 + eta)
n3 = (1.0/4.0) * (1.0 + eps) * (1.0 + eta)

x0, x1, x2, x3 = sym.symbols('x0 x1 x0 x1')
y0, y1, y2, y3 = sym.symbols('z0 z1 z0 z1')

interp = n0*x0 + n1*x1 + n2*x2 + n3*x3
objective_function = (x - (n0*x0 + n1*x1 + n2*x2 + n3*x3))

jacobian = sym.Matrix([sym.diff(objective_function, eps), sym.diff(objective_function, eta)])
hessian = sym.Matrix(([sym.diff(objective_function, eps, eps), sym.diff(objective_function, eps, eta)],
                      [sym.diff(objective_function, eta, eps), sym.diff(objective_function, eta, eta)]))

vertices = (4.9, 5.0, 4.9, 5.0)
solution = sym.solve(objective_function, eps)

sym.pprint(solution)
# sym.pprint(jacobian)
# print hessian
# print(hessian[0,1])

