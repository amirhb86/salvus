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

if __name__ == '__main__':
    main()
