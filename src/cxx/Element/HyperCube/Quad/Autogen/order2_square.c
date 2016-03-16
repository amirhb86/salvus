/******************************************************************************
 *                     Code generated with sympy 0.7.6.1                      *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                       This file is part of 'project'                       *
 ******************************************************************************/
#include "order2_square.h"
#include <math.h>

void interpolate_order2_square(double epsilon, double eta, double *out_5124881133433324070) {

   out_5124881133433324070[0] = 0.25*epsilon*eta*(epsilon - 1.0)*(eta - 1.0);
   out_5124881133433324070[1] = -0.5*eta*(epsilon - 1.0)*(epsilon + 1.0)*(eta - 1.0);
   out_5124881133433324070[2] = 0.25*epsilon*eta*(epsilon + 1.0)*(eta - 1.0);
   out_5124881133433324070[3] = -0.5*epsilon*(epsilon - 1.0)*(eta - 1.0)*(eta + 1.0);
   out_5124881133433324070[4] = 1.0*(epsilon - 1.0)*(epsilon + 1.0)*(eta - 1.0)*(eta + 1.0);
   out_5124881133433324070[5] = -0.5*epsilon*(epsilon + 1.0)*(eta - 1.0)*(eta + 1.0);
   out_5124881133433324070[6] = 0.25*epsilon*eta*(epsilon - 1.0)*(eta + 1.0);
   out_5124881133433324070[7] = -0.5*eta*(epsilon - 1.0)*(epsilon + 1.0)*(eta + 1.0);
   out_5124881133433324070[8] = 0.25*epsilon*eta*(epsilon + 1.0)*(eta + 1.0);

}

void interpolate_eps_derivative_order2_square(double epsilon, double eta, double *out_4655011469098782191) {

   out_4655011469098782191[0] = 0.25*epsilon*eta*(eta - 1.0) + 0.25*eta*(epsilon - 1.0)*(eta - 1.0);
   out_4655011469098782191[1] = -0.5*eta*(epsilon - 1.0)*(eta - 1.0) - 0.5*eta*(epsilon + 1.0)*(eta - 1.0);
   out_4655011469098782191[2] = 0.25*epsilon*eta*(eta - 1.0) + 0.25*eta*(epsilon + 1.0)*(eta - 1.0);
   out_4655011469098782191[3] = -0.5*epsilon*(eta - 1.0)*(eta + 1.0) - 0.5*(epsilon - 1.0)*(eta - 1.0)*(eta + 1.0);
   out_4655011469098782191[4] = 1.0*(epsilon - 1.0)*(eta - 1.0)*(eta + 1.0) + 1.0*(epsilon + 1.0)*(eta - 1.0)*(eta + 1.0);
   out_4655011469098782191[5] = -0.5*epsilon*(eta - 1.0)*(eta + 1.0) - 0.5*(epsilon + 1.0)*(eta - 1.0)*(eta + 1.0);
   out_4655011469098782191[6] = 0.25*epsilon*eta*(eta + 1.0) + 0.25*eta*(epsilon - 1.0)*(eta + 1.0);
   out_4655011469098782191[7] = -0.5*eta*(epsilon - 1.0)*(eta + 1.0) - 0.5*eta*(epsilon + 1.0)*(eta + 1.0);
   out_4655011469098782191[8] = 0.25*epsilon*eta*(eta + 1.0) + 0.25*eta*(epsilon + 1.0)*(eta + 1.0);

}

void interpolate_eta_derivative_order2_square(double epsilon, double eta, double *out_9086922723928993073) {

   out_9086922723928993073[0] = 0.25*epsilon*eta*(epsilon - 1.0) + 0.25*epsilon*(epsilon - 1.0)*(eta - 1.0);
   out_9086922723928993073[1] = -0.5*eta*(epsilon - 1.0)*(epsilon + 1.0) - 0.5*(epsilon - 1.0)*(epsilon + 1.0)*(eta - 1.0);
   out_9086922723928993073[2] = 0.25*epsilon*eta*(epsilon + 1.0) + 0.25*epsilon*(epsilon + 1.0)*(eta - 1.0);
   out_9086922723928993073[3] = -0.5*epsilon*(epsilon - 1.0)*(eta - 1.0) - 0.5*epsilon*(epsilon - 1.0)*(eta + 1.0);
   out_9086922723928993073[4] = 1.0*(epsilon - 1.0)*(epsilon + 1.0)*(eta - 1.0) + 1.0*(epsilon - 1.0)*(epsilon + 1.0)*(eta + 1.0);
   out_9086922723928993073[5] = -0.5*epsilon*(epsilon + 1.0)*(eta - 1.0) - 0.5*epsilon*(epsilon + 1.0)*(eta + 1.0);
   out_9086922723928993073[6] = 0.25*epsilon*eta*(epsilon - 1.0) + 0.25*epsilon*(epsilon - 1.0)*(eta + 1.0);
   out_9086922723928993073[7] = -0.5*eta*(epsilon - 1.0)*(epsilon + 1.0) - 0.5*(epsilon - 1.0)*(epsilon + 1.0)*(eta + 1.0);
   out_9086922723928993073[8] = 0.25*epsilon*eta*(epsilon + 1.0) + 0.25*epsilon*(epsilon + 1.0)*(eta + 1.0);

}

void closure_mapping_order2_square(int *out_5427770538133181134) {

   out_5427770538133181134[0] = 4;
   out_5427770538133181134[1] = 1;
   out_5427770538133181134[2] = 5;
   out_5427770538133181134[3] = 7;
   out_5427770538133181134[4] = 3;
   out_5427770538133181134[5] = 0;
   out_5427770538133181134[6] = 2;
   out_5427770538133181134[7] = 8;
   out_5427770538133181134[8] = 6;

}

void gll_weights_order2_square(double *out_1923134116974207604) {

   out_1923134116974207604[0] = 0.333333333333333;
   out_1923134116974207604[1] = 1.33333333333333;
   out_1923134116974207604[2] = 0.333333333333333;

}

void gll_coordinates_order2_square(double *out_2127667241337835057) {

   out_2127667241337835057[0] = -1.0;
   out_2127667241337835057[1] = 0.0;
   out_2127667241337835057[2] = 1.0;

}
