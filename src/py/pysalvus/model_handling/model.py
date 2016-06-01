# -*- coding: utf-8 -*-

import os
import shutil
import exodus

from sets import Set

import numpy as np

DEFAULT_BLOCK = 1
TIMESTEP_ZERO = 1


def getElemsAndNodes(input_file):
    "Gets the number of elements and nodes per elements in `input_file`"
    exo_from = exodus.exodus(input_file,"r",array_type='ctype')
    BLOCK_ID = 1
    _, _, nNodeElm, _ = exo_from.elem_blk_info(BLOCK_ID)
    num_elems = exo_from.num_elems()
    return (num_elems,nNodeElm)

def getExodusCopy(exodus_file_existing,exodus_file_new,parameter_names):
    '''
    Gets a copy of an existing exodus file, for the purpose of adding or changing existing parameter values (VP,VS,etc)
    :param exodus_file_existing Existing file
    :param exodus_file_new Exodus cannot write to an existing file, so write to a new one.
    :param parameter_names The variable names being added to the file. We cannot add these after the fact, so they are added now at file creation.
    :return The newly created exodus file
    '''
    
    exo_from = exodus.exodus(exodus_file_existing,"r",array_type='ctype')
    num_vars_from = exo_from.get_element_variable_number()
    
    if num_vars_from == 0:
        addElementVariables = parameter_names
        exo = exodus.copyTransfer(exodus_file_existing,exodus_file_new)
        print("Adding {} to no existing variables".format(addElementVariables))
        exo.set_element_variable_number(len(addElementVariables))
        for (i,name) in enumerate(addElementVariables):                        
            exo.put_element_variable_name(name,i+1)
        
        return exo
    else:
        existing_names = exo_from.get_element_variable_names()
        set_existing_names = Set(existing_names)
        set_parameter_names = Set(parameter_names)
        set_new_names = set_parameter_names - set_existing_names
        addElementVariables = list(set_new_names)
        print("Adding {} to existing variables {}".format(addElementVariables,existing_names))
        return exodus.copyTransfer(exodus_file_existing,exodus_file_new,additionalElementVariables=addElementVariables)


def setMaterialParameter(exoObj,parameter_names, parameter_values):
    # type: (exoObj: Any ,parameter_names: [str], parameter_values: [[double]])
    '''
    Sets value to an element parameter. The parameter name must already be in the file.
    :param exoObj Exodus file object
    :param paramter_names Ex: ["VP_0","VP_1","VP_2","fluid"]
    :param paramter_values Ex: For 3-element mesh "[[4,4,4],[4,4,4],[4,4,4],[0,0,0]]" (fluid=0 is elastic)
    '''
    BLOCK_ID = 1
    if not (parameter_names[0] in exoObj.get_element_variable_names()):
        raise Exception("ERROR: paramter name {} not in exodus file (has: {}) -- add it using getExodusCopy(,,parameter_names=names)!".format(parameter_names[0],exoObj.get_element_variable_names()))
    
    for (i,(parameter_name,parameter_value)) in enumerate(zip(parameter_names,parameter_values)):        
        exoObj.put_element_variable_values(BLOCK_ID, parameter_name, TIMESTEP_ZERO, parameter_value)
        
