# -*- coding: utf-8 -*-

from setuptools import setup

setup(
        name="pysalvus",
        version="0.1",
        py_modules=["pysalvus"],
        install_requires=[
            'Click'
        ],
        entry_points="""
        [console_scripts]
        pysalvus_code=pysalvus.pysalvus:code_generation
        pysalvus_model=pysalvus.pysalvus:model_handling
        pysalvus_solver=pysalvus.pysalvus:solver_operation
        pysalvus_optimize=pysalvus.pysalvus:optimization_tools
        """
)
