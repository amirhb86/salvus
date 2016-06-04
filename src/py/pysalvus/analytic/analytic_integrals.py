from __future__ import print_function

import click
import sympy


@click.group()
def main():
    pass


@main.command()
@click.argument('degree', default=4)
def analytic_2D_integral(degree):

    x, y = sympy.symbols('x y')
    xrange, yrange = range(degree+1), range(degree) + [degree-1]
    poly = ' + '.join(['x ** %d * y ** %d' % (ox, oy) for ox, oy in zip(xrange, yrange)])
    int = sympy.integrate(sympy.integrate(poly, (x, -2, 1)), (y, -6, 1))
    sympy.pprint(poly)
    print (int)

@main.command()
@click.argument('degree', default=4)
def analytic_1D_integral(degree):

    x, y = sympy.symbols('x y')
    xrange, yrange = range(degree+1), range(degree) + [degree-1]
    poly = ' + '.join(['x ** %d * y ** %d' % (ox, oy) for ox, oy in zip(xrange, yrange)])
    int = sympy.integrate(poly, (x, -2, 1))
    print ('INT_X', int.subs(y, 1))
    int = sympy.integrate(poly, (y, -6, 1))
    print ('INT_Y', int.subs(x, 1))

@main.command()
def analytic_acoustic_to_elastic():

    rho, v, x = sympy.symbols('rho v x')
    expr = rho * v
    expr = expr.subs(((rho, 2600), (v, 1)))
    print (sympy.integrate(expr, (x, 0, 50000)))


if __name__ == '__main__':
    main()
