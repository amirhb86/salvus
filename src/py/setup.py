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
        pysalvus=pysalvus.pysalvus:cli
        """
)
