#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests features from the code generation module.
"""

import numpy as np

from .code_generator import \
    generate_closure_mapping

def test_generate_closure_mapping():

    # reference closure for n = 4
    ref_closure = [6, 7, 8, 11, 12, 13, 16, 17, 18, 1, 2, 3, 9, 14,
                   19, 23, 22, 21, 15, 10, 5, 0, 4, 24, 20]
    
    closure = generate_closure_mapping(4)
    np.testing.assert_array_equal(closure, ref_closure)
