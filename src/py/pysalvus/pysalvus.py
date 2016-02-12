#! -*- coding: utf-8 -*-
from __future__ import print_function

import click
import numpy as np

from code_generation import code_generator
from model_handling import model


@click.group()
def cli():
    """
    pysalvus
    """
    pass


@cli.group()
def code_generation():
    """
    Test.
    """
    pass


@code_generation.command()
@click.option('--dimension', type=click.Choice(['2', '3']), required=True,
              help='Dimension of element.')
@click.option('--polynomial_order', required=True,
              help='Lagrange polynomial order.')
def hyperCube_generate_gll_basis(dimension, polynomial_order):
    """
    Autogenerate the appropriate code for a tensorized gll basis on a 2D/3D
    hypercube element.

    :param dimension: Dimension of element.
    :param polynomial_order: Lagrange polynomial order.
    """
    if dimension == '2':
        code_generator.tensorized_basis_2D(int(polynomial_order))
    else:
        raise NotImplementedError


@cli.group()
def model_handling():
    """
    Test.
    """
    pass


@model_handling.command()
@click.option('--input_file', help='Exodus file to add parameter to',
              type=click.Path(readable=True), required=True)
@click.option('--name', help='Name of parameter.', default='Velocity')
@click.option('--value', help='Value of constant parameter.', default=4.0)
@click.option('--output_file', help='Exodus file to to', default=None)
def add_constant_material_parameter(input_file, name, value, output_file):
    '''
    Test.
    '''
    working_model = model.ExodusModel(exodus_file=input_file)
    working_model.readFromExodus()

    parameter = np.ones(working_model.number_of_nodes) * value
    working_model.addMaterialParameter(name.lower(), parameter)

    working_model.write(output_file)


@cli.group()
def solver_operation():
    """
    Test.
    """
    pass


@cli.group()
def optimization_tools():
    """
    Test.
    """
    pass


if __name__ == "__main__":
    cli()
