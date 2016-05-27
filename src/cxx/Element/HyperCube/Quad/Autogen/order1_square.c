/******************************************************************************
 *                     Code generated with sympy 0.7.6.1                      *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                       This file is part of 'project'                       *
 ******************************************************************************/
#include "Element/HyperCube/Autogen/order1_square.h"
#include <math.h>

void interpolate_order1_square(double epsilon, double eta, double *out_8090246698051949025) {

   out_8090246698051949025[0] = 0.25*(epsilon - 1.0)*(eta - 1.0);
   out_8090246698051949025[1] = -0.25*(epsilon + 1.0)*(eta - 1.0);
   out_8090246698051949025[2] = -0.25*(epsilon - 1.0)*(eta + 1.0);
   out_8090246698051949025[3] = 0.25*(epsilon + 1.0)*(eta + 1.0);

}

void interpolate_eps_derivative_order1_square(double eta, double *out_9073736670949263119) {

   out_9073736670949263119[0] = 0.25*eta - 0.25;
   out_9073736670949263119[1] = -0.25*eta + 0.25;
   out_9073736670949263119[2] = -0.25*eta - 0.25;
   out_9073736670949263119[3] = 0.25*eta + 0.25;

}

void interpolate_eta_derivative_order1_square(double epsilon, double *out_8473474700611321857) {

   out_8473474700611321857[0] = 0.25*epsilon - 0.25;
   out_8473474700611321857[1] = -0.25*epsilon - 0.25;
   out_8473474700611321857[2] = -0.25*epsilon + 0.25;
   out_8473474700611321857[3] = 0.25*epsilon + 0.25;

}

void closure_mapping_order1_square(int *out_2820443818077165713) {

   out_2820443818077165713[0] = 0;
   out_2820443818077165713[1] = 1;
   out_2820443818077165713[2] = 3;
   out_2820443818077165713[3] = 2;

}

void gll_weights_order1_square(double *out_3122137466506344549) {

   out_3122137466506344549[0] = 1.0;
   out_3122137466506344549[1] = 1.0;

}

void gll_coordinates_order1_square(double *out_3093496787365923720) {

   out_3093496787365923720[0] = -1.0;
   out_3093496787365923720[1] = 1.0;

}
