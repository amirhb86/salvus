import os
import sympy as sym
from sympy.physics.quantum import TensorProduct
from sympy.utilities.codegen import CCodeGen
from quadrature_points_weights import \
    gauss_lobatto_legendre_quadruature_points_weights


def generating_polynomial_lagrange(N, coordinate, gll_points):
    """Symbolically calculates the generating polynomial for a lagrange polynomial of a given order.

    Our HyperCube elements use N + 1 Lagrange polynomials of order N for the spatial discretisation. This function
    returns N + 1 generating polynomials, which can be used to calculate the value of a Lagrange polynomial at
    an arbitrary location. This is particularily useful, because the full discretisation on a D dimensional element
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
    gll_coordinates, _ = gauss_lobatto_legendre_quadruature_points_weights(order + 1)
    basis = basis.subs([(v, c) for v, c in zip(eps_gll, gll_coordinates)])
    basis = basis.subs([(v, c) for v, c in zip(eta_gll, gll_coordinates)])
    sym.pprint(basis)

    # Get gradient of basis functions.
    basis_gradient_eps = sym.Matrix([sym.diff(i, eps) for i in basis])
    basis_gradient_eta = sym.Matrix([sym.diff(i, eta) for i in basis])

    # Get diagonal mass matrix
    mass_matrix = rho * basis * basis.T
    mass_matrix_diagonal = sym.Matrix([mass_matrix[i,i] for i in range(total_integration_points)])

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
        'diagonal_mass_matrix_order{}_square'.format(order), mass_matrix_diagonal,
        argument_sequence=None))
    autocode.write(routines, 'order{}_square'.format(order), to_files=True)
