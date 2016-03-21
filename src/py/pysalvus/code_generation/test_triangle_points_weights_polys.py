#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests features from the code generation module.
"""

import numpy as np
import sympy as sym

from .triangle_points_weights_polys import p3ElementPolynomials,p3ReferenceNodes
    
def test_p3_triangle_phi():

    r,s = sym.symbols("r,s")
    rsn = p3ReferenceNodes()
    Np = len(rsn)
    phi = p3ElementPolynomials(rsn,r,s)

    for (i,phi_i) in enumerate(phi):
        
        # test
        phi_if = sym.lambdify((r,s),phi_i)
        b = np.zeros(Np)
        b[i] = 1.0

        # Check if 1 at *this* point rsn[i], and 0 at other points
        # fulfilling lagrange polynomial promise
        lagrange_eval = np.asarray([phi_if(rn,sn) for (rn,sn) in rsn])
        np.testing.assert_allclose(lagrange_eval,b,atol=1e-13)

        
