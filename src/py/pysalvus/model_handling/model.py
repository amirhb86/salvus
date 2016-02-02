# -*- coding: utf-8 -*-

import os
import shutil
import exodus

TIMESTEP_ZERO = 1
TEMP_EXODUS_FILENAME = '.tmp.exo2'

class ExodusModel(object):
    def __init__(self, exodus_file):
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
        self.exodus_template = exodus.exodus(self.exodus_file, array_type='numpy')

    def __reload_from_tmp(self):
        self.exodus_template = exodus.exodus(TEMP_EXODUS_FILENAME, array_type='numpy')

    def readFromExodus(self):
        self.__read()

    def write(self, file_name):
        if file_name:
            tmp = self.exodus_template.copy(file_name)
            tmp.close()
        else:
            os.remove(self.exodus_file)
            shutil.copy(TEMP_EXODUS_FILENAME, self.exodus_file)


    def addMaterialParameter(self, parameter_name, parameter_value):
        if parameter_name in self.material_parameters:
            raise RuntimeError("Parameter already exists. Use updateMaterialParameter() to overwrite.")
        else:
            buffer = self.exodus_template.copy(self.temporary_exodus_file)
            buffer.set_node_variable_number(self.number_of_nodal_variables + 1)
            buffer.put_node_variable_name(parameter_name, self.number_of_nodal_variables + 1)
            buffer.put_node_variable_values(parameter_name, TIMESTEP_ZERO, parameter_value)
            buffer.close()
            self.__reload_from_tmp()