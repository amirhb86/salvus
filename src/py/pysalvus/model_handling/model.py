# -*- coding: utf-8 -*-

import os
import shutil
import exodus

DEFAULT_BLOCK = 1
TIMESTEP_ZERO = 1

def getExodusCopy(exodus_file_existing,exodus_file_new):
    return exodus.copyTransfer(exodus_file_existing,exodus_file_new)
    
def addMaterialParameter(exoObj,parameter_name, parameter_values):
    BLOCK_ID = 1
    if parameter_name in exoObj.get_node_variable_names():        
        raise RuntimeError("Parameter already exists. Use updateMaterialParameter() to overwrite.")
    else:
        _, _, nNodeElm, _ = exoObj.elem_blk_info(BLOCK_ID)
        exoObj.set_element_variable_number(nNodeElm)
        for node in range(nNodeElm):
            var_name = "{}_{}".format(parameter_name, node)
            exoObj.put_element_variable_name(var_name, node+1)
            exoObj.put_element_variable_values(BLOCK_ID, var_name, TIMESTEP_ZERO, parameter_values)

def updateMaterialParameter(exoObj,parameter_name, parameter_values):
    BLOCK_ID = DEFAULT_BLOCK
    if "{}_0".format(parameter_name) in exoObj.get_element_variable_names():        
        _, _, nNodeElm, _ = exoObj.elem_blk_info(BLOCK_ID)
        if nNodeElm != exoObj.get_element_variable_number():
            raise RuntimeError("Number of element parameters and nodes doesn't match!")
        for node in range(nNodeElm):
            var_name = "{}_{}".format(parameter_name, node)
            exoObj.put_element_variable_values(BLOCK_ID, var_name, TIMESTEP_ZERO, parameter_values)        
    else:
        raise RuntimeError("Parameter doesn't exist. Use addMaterialParameter() to add it.")
            
