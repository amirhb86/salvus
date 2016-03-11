import os

import io
import sympy as sym
from sympy.physics.quantum import TensorProduct
from sympy.utilities.codegen import CCodeGen
from quadrature_points_weights import \
    gauss_lobatto_legendre_quadruature_points_weights


def generating_polynomial_lagrange(N, coordinate, gll_points):
    """Symbolically calculates the generating polynomial for a lagrange polynomial of a given order.

    Our HyperCube elements use N + 1 Lagrange polynomials of order N for the spatial discretisation. This function
    returns N + 1 generating polynomials, which can be used to calculate the value of a Lagrange polynomial at
    an arbitrary location. This is particularly useful, because the full discretisation on a D dimensional element
    can be obtained by the tensor product these polynomials.

    This routine is derived on pg. 304 of Andreas' book.

    :param N: Order of Lagrange polynomial.
    :param coordinate: Coordinate in reference element (i.e. eps, eta, phi).
    :returns: A N + 1 element sympy vector of generating polynomials of order N.
    """

    # Symbols.
    coord, phi, d_phi = sym.symbols('{} Phi dPhi'.format(coordinate))

    # Equation A.19
    phi = 1
    for x in gll_points:
        phi *= (coord - x)

    # Get derivative in coordinate direction
    d_phi = sym.diff(phi, coord)

    # Equation A.20
    generators = []
    for x in gll_points:
        generators.append((1 / (coord - x)) * (phi / d_phi.subs(coord, x)))

    return sym.Matrix(generators)

def generate_closure_mapping(N):
    """Generates the mapping from PETSc's element closure to one we definde.

    The DMPlexVecGet(/Set)Closure assumes a certain ordering of the element unknowns. We want to reorder
    these unknowns into some arbitrary order, perhaps one that may be convienient for use in a tensor
    product space. This function generates the mapping to a 2-D tensorized basis as below:
    
                       (edge_2)
        v(n*(n+1))_________________ v((n+1)*(n+1)-1)
                  |               |
                 .|               | .
                 .|               | .               (edge_1)
 (edge_3)        .|               | .
        v(2*(n+1))|               | v(3*(n+1)-1)
                  |               |
            v(n+1)|_______________| v(2*(n+1)-1)

                  v0, v1, ...,   vn

                       (edge_0)


        (j)
        ^
        |

        |_ _ _ > (i)

    :param N: Order of Lagrange basis.
    """

    # Dofs per component.
    numPerVertex = 1
    numPerEdge = N - 1
    numPerFace = (N - 1) ** 2

    vertices = [0, N, (N+1)*(N+1)-1, N*(N+1)]
    face = [i for j in range(1,N) for i in range(j*(N+1)+1,j*(N+1)+N)]
    edge_0 = range(1, N)
    edge_1 = range(2*(N+1)-1, (N+1)*(N+1)-1, N+1)
    edge_2 = range((N+1)*(N+1)-2, N*(N+1), -1)
    edge_3 = range((N-1)*(N+1), 0, -1*(N+1))

    return face + edge_0 + edge_1 + edge_2 + edge_3 + vertices

def tensorized_basis_2D(order):

    total_integration_points = (order + 1) * (order + 1)
    eps, eta, rho = sym.symbols('epsilon eta rho')
    eps_gll = sym.symbols('epsilon_0:%d' % (order + 1))
    eta_gll = sym.symbols('eta_0:%d' % (order + 1))

    # Get N + 1 lagrange polynomials in each direction.
    generator_eps = generating_polynomial_lagrange(order, 'epsilon', eps_gll)
    generator_eta = generating_polynomial_lagrange(order, 'eta', eta_gll)

    # Get tensorized basis.
    basis = TensorProduct(generator_eta, generator_eps)
    gll_coordinates, gll_weights = gauss_lobatto_legendre_quadruature_points_weights(order + 1)
    basis = basis.subs([(v, c) for v, c in zip(eps_gll, gll_coordinates)])
    basis = basis.subs([(v, c) for v, c in zip(eta_gll, gll_coordinates)])

    # Get gradient of basis functions.
    basis_gradient_eps = sym.Matrix([sym.diff(i, eps) for i in basis])
    basis_gradient_eta = sym.Matrix([sym.diff(i, eta) for i in basis])

    # Get closure mapping.
    closure = sym.Matrix(generate_closure_mapping(order), dtype=int)

    # Write code
    routines = []
    autocode = CCodeGen()
    routines.append(autocode.routine(
        'interpolate_order{}_square'.format(order), basis,
        argument_sequence=None))
    routines.append(autocode.routine(
        'interpolate_eps_derivative_order{}_square'.format(order), basis_gradient_eps,
        argument_sequence=None))
    routines.append(autocode.routine(
        'interpolate_eta_derivative_order{}_square'.format(order), basis_gradient_eta,
        argument_sequence=None))
    routines.append(autocode.routine(
        'closure_mapping_order{}_square'.format(order), closure,
        argument_sequence=None))
    routines.append(autocode.routine(
        'gll_weights_order{}_square'.format(order), sym.Matrix(gll_weights),
        argument_sequence=None))
    routines.append(autocode.routine(
        'gll_coordinates_order{}_square'.format(order), sym.Matrix(gll_coordinates),
        argument_sequence=None))
    autocode.write(routines, 'order{}_square'.format(order), to_files=True)

    # reformat some code.
    for code, lend in zip(['order{}_square.c', 'order{}_square.h'], [' {', ';']):
        with io.open(code.format(order), 'rt') as fh:
            text = fh.readlines()
            text = [line.replace('double', 'int') if 'closure' in line else line for line in text]

        with io.open(code.format(order), 'wt') as fh:
            fh.writelines(text)

