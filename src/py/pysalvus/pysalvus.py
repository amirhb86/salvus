#! -*- coding: utf-8 -*-
from __future__ import print_function
from model_handling import exodus as elib

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
        code_generator.tensorized_basis_3D(int(polynomial_order))


@cli.group()
def model_handling():
    """
    Test.
    """
    pass


# minimum elastic parameters: RHO, C11, C33, C55, C13
# Usage: python pysalvus.py model_handling add_constant_material_parameter --input_file ../../../../salvus_data/test_meshes_tri/simple_quadmesh_2x2.e --name VP --value 4 --output_file simple_quadmesh_2x2.param2.e
# Usage: python pysalvus.py model_handling add_constant_material_parameter --input-file simple_quadmesh_2x2.e --name RHO --value 1.0 --name C11 --value 4.0 --name C33 --value 4.0 --name C55 --value 4.0 --name C13 --value 4.0 --output-file simple_quadmesh_2x2.C_all.e
# 
@model_handling.command()
@click.option('--input-file', help='Exodus file to add parameter to',
              type=click.Path(readable=True), required=True)
@click.option('--name', help='Name of parameter.', default='VP',multiple=True)
@click.option('--value', help='Value of constant parameter.', default=4.0,multiple=True)
@click.option('--output-file', help='Exodus file to to', required=True)
def add_constant_material_parameter(input_file, name, value, output_file):
    '''
    Adds single parameter with a constant value for entire mesh.
    > python pysalvus.py model_handling add_constant_material_parameter --input_file mesh-in.e --name VP --value 4 --output_file mesh-out.e
    '''    
    (num_elems,nNodeElm) = model.getElemsAndNodes(input_file)
    parameter_values = []
    parameter_names = []
    for (this_name,this_value) in zip(name,value):        
        for i in range(0,nNodeElm):
            parameter_names.append("{}_{}".format(this_name,i))
            parameter_values.append(np.ones(num_elems) * this_value)
                    
    exo = model.getExodusCopy(input_file,output_file,parameter_names)
    model.setMaterialParameter(exo,parameter_names,parameter_values)
    exo.close()

# Usage: python pysalvus.py model_handling add_fluid_flag --input-file simple_quadmesh_2x2.e --fluid-type "elastic" --output-file simple_quadmesh_2x2.elastic.e    
@model_handling.command()
@click.option('--input-file', help='Exodus file to add parameter to',
             type=click.Path(readable=True), required=True)
@click.option('--fluid-type', help='fluid or elastic solid?.', default="fluid")
@click.option('--output-file', help='Exodus file to add parameter to',
              required=True)
def add_fluid_flag(input_file, fluid_type, output_file):
    '''
    Just a quick function to mark everything as fluid.
    :param input_file:
    :param output_file:
    :return:
    '''
    (num_elems,nNodeElm) = model.getElemsAndNodes(input_file)    
    
    if fluid_type == "fluid":
        values = np.ones(num_elems)
    elif fluid_type == "elastic":
        values = np.zeros(num_elems)

    exo = model.getExodusCopy(input_file,output_file,["fluid"])
    model.setMaterialParameter(exo,["fluid"],[values])
    exo.close()


    
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
