/******************************************************************************
 *                     Code generated with sympy 0.7.6.1                      *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                       This file is part of 'project'                       *
 ******************************************************************************/
#include "order3_square.h"
#include <math.h>

void interpolate_order3_square(double epsilon, double eta, double *out_5520632088669242493) {

   out_5520632088669242493[0] = 0.390625*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958);
   out_5520632088669242493[1] = -0.873464053710856*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958);
   out_5520632088669242493[2] = 0.873464053710855*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958);
   out_5520632088669242493[3] = -0.390625*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958);
   out_5520632088669242493[4] = -0.873464053710855*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0);
   out_5520632088669242493[5] = 1.953125*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0);
   out_5520632088669242493[6] = -1.953125*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0);
   out_5520632088669242493[7] = 0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0);
   out_5520632088669242493[8] = 0.873464053710855*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_5520632088669242493[9] = -1.953125*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_5520632088669242493[10] = 1.953125*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_5520632088669242493[11] = -0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_5520632088669242493[12] = -0.390625*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0);
   out_5520632088669242493[13] = 0.873464053710856*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0);
   out_5520632088669242493[14] = -0.873464053710855*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0);
   out_5520632088669242493[15] = 0.390625*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0);

}

void interpolate_eps_derivative_order3_square(double epsilon, double eta, double *out_1735109083759235119) {

   out_1735109083759235119[0] = 0.390625*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958) + 0.390625*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958) + 0.390625*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958);
   out_1735109083759235119[1] = -0.873464053710856*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958) - 0.873464053710856*(epsilon - 1.0)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958) - 0.873464053710856*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958);
   out_1735109083759235119[2] = 0.873464053710855*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958) + 0.873464053710855*(epsilon - 1.0)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958) + 0.873464053710855*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958);
   out_1735109083759235119[3] = -0.390625*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958) - 0.390625*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958) - 0.390625*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958);
   out_1735109083759235119[4] = -0.873464053710855*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0) - 0.873464053710855*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0) - 0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0);
   out_1735109083759235119[5] = 1.953125*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0) + 1.953125*(epsilon - 1.0)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0) + 1.953125*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0);
   out_1735109083759235119[6] = -1.953125*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0) - 1.953125*(epsilon - 1.0)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0) - 1.953125*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0);
   out_1735109083759235119[7] = 0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0) + 0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0) + 0.873464053710855*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958)*(eta + 1.0);
   out_1735109083759235119[8] = 0.873464053710855*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0) + 0.873464053710855*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0) + 0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_1735109083759235119[9] = -1.953125*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0) - 1.953125*(epsilon - 1.0)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0) - 1.953125*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_1735109083759235119[10] = 1.953125*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0) + 1.953125*(epsilon - 1.0)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0) + 1.953125*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_1735109083759235119[11] = -0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0) - 0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0) - 0.873464053710855*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_1735109083759235119[12] = -0.390625*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0) - 0.390625*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0) - 0.390625*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0);
   out_1735109083759235119[13] = 0.873464053710856*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0) + 0.873464053710856*(epsilon - 1.0)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0) + 0.873464053710856*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0);
   out_1735109083759235119[14] = -0.873464053710855*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0) - 0.873464053710855*(epsilon - 1.0)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0) - 0.873464053710855*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0);
   out_1735109083759235119[15] = 0.390625*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0) + 0.390625*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0) + 0.390625*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0);

}

void interpolate_eta_derivative_order3_square(double epsilon, double eta, double *out_7408361074971146125) {

   out_7408361074971146125[0] = 0.390625*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958) + 0.390625*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta + 0.447213595499958) + 0.390625*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 0.447213595499958)*(eta + 0.447213595499958);
   out_7408361074971146125[1] = -0.873464053710856*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958) - 0.873464053710856*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958) - 0.873464053710856*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958);
   out_7408361074971146125[2] = 0.873464053710855*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958) + 0.873464053710855*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958) + 0.873464053710855*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958);
   out_7408361074971146125[3] = -0.390625*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958) - 0.390625*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958) - 0.390625*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958);
   out_7408361074971146125[4] = -0.873464053710855*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta - 0.447213595499958) - 0.873464053710855*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta + 1.0) - 0.873464053710855*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 0.447213595499958)*(eta + 1.0);
   out_7408361074971146125[5] = 1.953125*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958) + 1.953125*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 1.0) + 1.953125*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 1.0);
   out_7408361074971146125[6] = -1.953125*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958) - 1.953125*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 1.0) - 1.953125*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 1.0);
   out_7408361074971146125[7] = 0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta - 0.447213595499958) + 0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 1.0) + 0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 1.0);
   out_7408361074971146125[8] = 0.873464053710855*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta + 0.447213595499958) + 0.873464053710855*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 1.0)*(eta + 1.0) + 0.873464053710855*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0);
   out_7408361074971146125[9] = -1.953125*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958) - 1.953125*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 1.0) - 1.953125*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_7408361074971146125[10] = 1.953125*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958) + 1.953125*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 1.0) + 1.953125*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_7408361074971146125[11] = -0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 0.447213595499958) - 0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 1.0)*(eta + 1.0) - 0.873464053710855*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_7408361074971146125[12] = -0.390625*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 0.447213595499958)*(eta + 0.447213595499958) - 0.390625*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta - 0.447213595499958)*(eta + 1.0) - 0.390625*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(eta + 0.447213595499958)*(eta + 1.0);
   out_7408361074971146125[13] = 0.873464053710856*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958) + 0.873464053710856*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 1.0) + 0.873464053710856*(epsilon - 1.0)*(epsilon - 0.447213595499958)*(epsilon + 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_7408361074971146125[14] = -0.873464053710855*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958) - 0.873464053710855*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 1.0) - 0.873464053710855*(epsilon - 1.0)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta + 0.447213595499958)*(eta + 1.0);
   out_7408361074971146125[15] = 0.390625*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 0.447213595499958) + 0.390625*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta - 0.447213595499958)*(eta + 1.0) + 0.390625*(epsilon - 0.447213595499958)*(epsilon + 0.447213595499958)*(epsilon + 1.0)*(eta + 0.447213595499958)*(eta + 1.0);

}

void diagonal_mass_matrix_order3_square(double epsilon, double eta, double rho, double *out_8631827086440990231) {

   out_8631827086440990231[0] = 0.152587890625*rho*pow(epsilon - 1.0, 2)*pow(epsilon - 0.447213595499958, 2)*pow(epsilon + 0.447213595499958, 2)*pow(eta - 1.0, 2)*pow(eta - 0.447213595499958, 2)*pow(eta + 0.447213595499958, 2);
   out_8631827086440990231[1] = 0.762939453125*rho*pow(epsilon - 1.0, 2)*pow(epsilon - 0.447213595499958, 2)*pow(epsilon + 1.0, 2)*pow(eta - 1.0, 2)*pow(eta - 0.447213595499958, 2)*pow(eta + 0.447213595499958, 2);
   out_8631827086440990231[2] = 0.762939453125*rho*pow(epsilon - 1.0, 2)*pow(epsilon + 0.447213595499958, 2)*pow(epsilon + 1.0, 2)*pow(eta - 1.0, 2)*pow(eta - 0.447213595499958, 2)*pow(eta + 0.447213595499958, 2);
   out_8631827086440990231[3] = 0.152587890625*rho*pow(epsilon - 0.447213595499958, 2)*pow(epsilon + 0.447213595499958, 2)*pow(epsilon + 1.0, 2)*pow(eta - 1.0, 2)*pow(eta - 0.447213595499958, 2)*pow(eta + 0.447213595499958, 2);
   out_8631827086440990231[4] = 0.762939453125*rho*pow(epsilon - 1.0, 2)*pow(epsilon - 0.447213595499958, 2)*pow(epsilon + 0.447213595499958, 2)*pow(eta - 1.0, 2)*pow(eta - 0.447213595499958, 2)*pow(eta + 1.0, 2);
   out_8631827086440990231[5] = 3.814697265625*rho*pow(epsilon - 1.0, 2)*pow(epsilon - 0.447213595499958, 2)*pow(epsilon + 1.0, 2)*pow(eta - 1.0, 2)*pow(eta - 0.447213595499958, 2)*pow(eta + 1.0, 2);
   out_8631827086440990231[6] = 3.814697265625*rho*pow(epsilon - 1.0, 2)*pow(epsilon + 0.447213595499958, 2)*pow(epsilon + 1.0, 2)*pow(eta - 1.0, 2)*pow(eta - 0.447213595499958, 2)*pow(eta + 1.0, 2);
   out_8631827086440990231[7] = 0.762939453125*rho*pow(epsilon - 0.447213595499958, 2)*pow(epsilon + 0.447213595499958, 2)*pow(epsilon + 1.0, 2)*pow(eta - 1.0, 2)*pow(eta - 0.447213595499958, 2)*pow(eta + 1.0, 2);
   out_8631827086440990231[8] = 0.762939453125*rho*pow(epsilon - 1.0, 2)*pow(epsilon - 0.447213595499958, 2)*pow(epsilon + 0.447213595499958, 2)*pow(eta - 1.0, 2)*pow(eta + 0.447213595499958, 2)*pow(eta + 1.0, 2);
   out_8631827086440990231[9] = 3.814697265625*rho*pow(epsilon - 1.0, 2)*pow(epsilon - 0.447213595499958, 2)*pow(epsilon + 1.0, 2)*pow(eta - 1.0, 2)*pow(eta + 0.447213595499958, 2)*pow(eta + 1.0, 2);
   out_8631827086440990231[10] = 3.814697265625*rho*pow(epsilon - 1.0, 2)*pow(epsilon + 0.447213595499958, 2)*pow(epsilon + 1.0, 2)*pow(eta - 1.0, 2)*pow(eta + 0.447213595499958, 2)*pow(eta + 1.0, 2);
   out_8631827086440990231[11] = 0.762939453125*rho*pow(epsilon - 0.447213595499958, 2)*pow(epsilon + 0.447213595499958, 2)*pow(epsilon + 1.0, 2)*pow(eta - 1.0, 2)*pow(eta + 0.447213595499958, 2)*pow(eta + 1.0, 2);
   out_8631827086440990231[12] = 0.152587890625*rho*pow(epsilon - 1.0, 2)*pow(epsilon - 0.447213595499958, 2)*pow(epsilon + 0.447213595499958, 2)*pow(eta - 0.447213595499958, 2)*pow(eta + 0.447213595499958, 2)*pow(eta + 1.0, 2);
   out_8631827086440990231[13] = 0.762939453125*rho*pow(epsilon - 1.0, 2)*pow(epsilon - 0.447213595499958, 2)*pow(epsilon + 1.0, 2)*pow(eta - 0.447213595499958, 2)*pow(eta + 0.447213595499958, 2)*pow(eta + 1.0, 2);
   out_8631827086440990231[14] = 0.762939453125*rho*pow(epsilon - 1.0, 2)*pow(epsilon + 0.447213595499958, 2)*pow(epsilon + 1.0, 2)*pow(eta - 0.447213595499958, 2)*pow(eta + 0.447213595499958, 2)*pow(eta + 1.0, 2);
   out_8631827086440990231[15] = 0.152587890625*rho*pow(epsilon - 0.447213595499958, 2)*pow(epsilon + 0.447213595499958, 2)*pow(epsilon + 1.0, 2)*pow(eta - 0.447213595499958, 2)*pow(eta + 0.447213595499958, 2)*pow(eta + 1.0, 2);

}