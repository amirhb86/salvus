from __future__ import print_function

import click
import sympy
from ..code_generation.code_generator import generating_polynomial_lagrange
from ..code_generation.quadrature_points_weights import \
    gauss_lobatto_legendre_quadruature_points_weights


def poly_2D(degree, verbose=False):
    """Generate a 2D polynomial for analytical testing."""

    r, s = sympy.symbols('r s')
    xrange, yrange = range(degree + 1), range(degree) + [degree - 1]
    poly = ' + '.join(['r ** %d * s ** %d' % (ox, oy) for ox, oy in zip(xrange, yrange)])
    int = sympy.integrate(sympy.integrate(poly, (r, -2, 1)), (s, -6, 1))
    return poly


@main.command()
@click.argument('degree', default=4)
def analytic_2D_integral(degree):
    """Print the analytical polynomial for a given degree."""

    print("Polynomial: ")
    sympy.pprint(poly_2D(degree))


@main.command()
@click.argument('degree', default=4)
def gradient_of_2d_function(degree):
    """Print the gradient the testing polynomial.

    For a given degree, print the gradient of a 2D polynomial in the
    reference (r,s) direction.
    """
    r, s = sympy.symbols('r s')
    poly = poly_2D(degree)
    print("Polynomial: ")
    sympy.pprint(poly)
    print("Derivative r: ")
    sympy.pprint(sympy.diff(poly, x))
    print("Derivative s: ")
    sympy.pprint(sympy.diff(poly, y))


@main.command()
@click.argument('degree', default=4)
def grad_test_times_grad_field(degree):
    """Calculate \int \nabla \phi _i \cdot \nabla \phi _j."""

    # Setup variables.
    r, s = sympy.symbols('r s')
    r_gll = sympy.symbols('r_0:%d' % (degree + 1))
    s_gll = sympy.symbols('s_0:%d' % (degree + 1))
    total_integration_points = (degree + 1) * (degree + 1)

    # Get N + 1 lagrange polynomials in each direction.
    generator_r = generating_polynomial_lagrange(degree, 'r', r_gll)
    generator_s = generating_polynomial_lagrange(degree, 's', s_gll)

    # Evalute bases
    basis = sympy.physics.quantum.TensorProduct(generator_r, generator_s)
    gll_coordinates, gll_weights = gauss_lobatto_legendre_quadruature_points_weights(degree + 1)
    basis = basis.subs([(v, c) for v, c in zip(r_gll, gll_coordinates)])
    basis = basis.subs([(v, c) for v, c in zip(s_gll, gll_coordinates)])

    # Get gradient of basis functions.
    basis_gradient_r = [sympy.diff(i, r) for i in basis]
    basis_gradient_s = [sympy.diff(i, s) for i in basis]

    # Get interpolating field.
    field = poly_2D(degree)

    # Take gradient of this field.
    grad_field = sympy.Matrix([[sympy.diff(field, r).subs({r: i, s: j})
                                for j in gll_coordinates
                                for i in gll_coordinates],
                               [sympy.diff(field, s).subs({r: i, s: j})
                                for j in gll_coordinates
                                for i in gll_coordinates]])

    # Take gradient of the tensor basis.
    grad_test = sympy.Matrix([basis_gradient_r,
                              basis_gradient_s])

    # Compute \nabla \phi _i \cdot \nabla \phi _j
    grad_test_grad_field = [grad_test.col(i).dot(grad_field.col(j))
                            for i in range(total_integration_points)
                            for j in range(total_integration_points)]

    # Sum over \phi _i
    final_grad = [sum(grad_test_grad_field[i:i + total_integration_points])
                  for i in range(0, len(grad_test_grad_field), total_integration_points)]

    # Integrate over 2D domain.
    integrated = [sympy.integrate(f, (r, -1, 1))
                  for f in final_grad]
    integrated = [sympy.integrate(f, (s, -1, 1))
                  for f in integrated]

    print("NABLA PHI_i \cdot NABLA PHI_j")
    sympy.pprint(integrated)
    print(sum(integrated))


@main.command()
def analytic_acoustic_to_elastic():
    rho, v, x = sympy.symbols('rho v x')
    expr = rho * v
    expr = expr.subs(((rho, 2600), (v, 1)))
    print(sympy.integrate(expr, (x, 0, 50000)))


@click.group()
def main():
    pass


if __name__ == '__main__':
    main()
