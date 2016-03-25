# -*- coding: utf-8 -*-

import os
import shutil
import exodus

DEFAULT_BLOCK = 1
TIMESTEP_ZERO = 1
TEMP_EXODUS_FILENAME = '.tmp.exo2'

class ExodusModel(object):
    def __init__(self, exodus_file):
        self.additionalNodalParameter = []
        self.additionalNodalParameterValues = []
        self.exodus_file = exodus_file

    def __repr__(self):
        return "ExodusModel({})".format(self.exodus_file)

    def __del__(self):
        if os.path.exists(TEMP_EXODUS_FILENAME):
            os.remove(TEMP_EXODUS_FILENAME)
        self.exodus_template.close()

    @property
    def number_of_elements(self):
        return self.exodus_template.num_elems()

    @property
    def number_of_nodes(self):
        return self.exodus_template.num_nodes()

    @property
    def material_parameters(self):
        return self.exodus_template.get_node_variable_names()

    @property
    def number_of_nodal_variables(self):
        return self.exodus_template.get_node_variable_number()

    @property
    def temporary_exodus_file(self):
        if os.path.exists(TEMP_EXODUS_FILENAME):
            os.remove(TEMP_EXODUS_FILENAME)
        return TEMP_EXODUS_FILENAME

    def __read(self):
        exodus.copyTransfer(self.exodus_file, self.temporary_exodus_file)
        self.__reload_from_tmp()

    def __reload_from_tmp(self):
        self.exodus_template = exodus.exodus(TEMP_EXODUS_FILENAME, array_type='numpy', mode='r')

    def readFromExodus(self):
        self.__read()

    def write(self, file_name):
        if file_name:
            if os.path.exists(file_name):
                os.remove(file_name)
            tmp = exodus.copyTransfer(TEMP_EXODUS_FILENAME, file_name,
                                      additionalNodalVariables=self.additionalNodalParameter)
            _, _, nNodeElm, _ = tmp.elem_blk_info(DEFAULT_BLOCK)
            tmp.set_element_variable_number(nNodeElm)
            for name, value in zip(self.additionalNodalParameter, self.additionalNodalParameterValues):
                for node in range(nNodeElm):
                    new_evar_index = node + 1
                    var_name = "{}_{}".format(name, node)
                    tmp.put_element_variable_name(var_name, new_evar_index)
                    tmp.put_element_variable_values(DEFAULT_BLOCK, var_name, TIMESTEP_ZERO, value)
            tmp.close()
        else:
            os.remove(self.exodus_file)
            shutil.copy(TEMP_EXODUS_FILENAME, self.exodus_file)


    def addMaterialParameter(self, parameter_name, parameter_value):
        if parameter_name in self.material_parameters:
            raise RuntimeError("Parameter already exists. Use updateMaterialParameter() to overwrite.")
        else:
            self.additionalNodalParameter.append(parameter_name)
            self.additionalNodalParameterValues.append(parameter_value)
