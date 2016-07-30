#include "catch.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <Eigen/Dense>
#include <petsc.h>

#include <Problem/Problem.h>
#include <Mesh/Mesh.h>
#include <Model/ExodusModel.h>

#include <Element/Element.h>
#include <Element/ElementAdapter.h>
#include <Element/HyperCube/QuadP1.h>
#include <Element/HyperCube/TensorQuad.h>

#include <Utilities/Types.h>
#include <Utilities/Options.h>

extern "C" {
#include <Element/HyperCube/Autogen/quad_autogen.h>
}

using namespace std;
using namespace Eigen;

RealMat derivative4order(const PetscReal r, const PetscReal s, const PetscInt order) {
  PetscInt size = (order + 1) * (order + 1);
  RealVec test_r(size), test_s(size);
  if (order == 1) {
    interpolate_eps_derivative_order1_square(s,    test_r.data());
    interpolate_eta_derivative_order1_square(r,    test_s.data());
  } else if (order == 2) {
    interpolate_eps_derivative_order2_square(r, s, test_r.data());
    interpolate_eta_derivative_order2_square(r, s, test_s.data());
  } else if (order == 3) {
    interpolate_eps_derivative_order3_square(r, s, test_r.data());
    interpolate_eta_derivative_order3_square(r, s, test_s.data());
  } else if (order == 4) {
    interpolate_eps_derivative_order4_square(r, s, test_r.data());
    interpolate_eta_derivative_order4_square(r, s, test_s.data());
  } else if (order == 5) {
    interpolate_eps_derivative_order5_square(r, s, test_r.data());
    interpolate_eta_derivative_order5_square(r, s, test_s.data());
  } else if (order == 6) {
    interpolate_eps_derivative_order6_square(r, s, test_r.data());
    interpolate_eta_derivative_order6_square(r, s, test_s.data());
  } else if (order == 7) {
    interpolate_eps_derivative_order7_square(r, s, test_r.data());
    interpolate_eta_derivative_order7_square(r, s, test_s.data());
  } else if (order == 8) {
    interpolate_eps_derivative_order8_square(r, s, test_r.data());
    interpolate_eta_derivative_order8_square(r, s, test_s.data());
  } else if (order == 9) {
    interpolate_eps_derivative_order9_square(r, s, test_r.data());
    interpolate_eta_derivative_order9_square(r, s, test_s.data());
  } else if (order == 10) {
    interpolate_eps_derivative_order10_square(r, s, test_r.data());
    interpolate_eta_derivative_order10_square(r, s, test_s.data());
  }
  RealMat ret(size, 2);
  ret.col(0) = test_r; ret.col(1) = test_s;
  return ret;
}

TEST_CASE("Test tensor quad", "[tensor_quad]") {

  std::string e_file = "fluid_layer_over_elastic_cartesian_2D_50s.e";

  PetscOptionsClear(NULL);
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--mesh-file", e_file.c_str(),
      "--model-file", e_file.c_str(),
      "--polynomial-order", "4",
      NULL
  };

  // Fake setting via command line.
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(NULL, &argc, &argv, NULL);

  // Default vertices.
  QuadVtx vtx;
  vtx << -1, -1, +1, -1, +1, +1, -1, +1;

  // Initialize dummy mesh and model.
  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  std::unique_ptr<Problem>      problem(Problem::Factory(options));
  std::unique_ptr<ExodusModel>  model(new ExodusModel(options));
  std::unique_ptr<Mesh>         mesh(new Mesh(options));

  mesh->read();
  model->read();
  mesh->setupGlobalDof(model, options);

  /* Loop over all polynomial orders. */
  for (PetscInt i = 1; i < TensorQuad<QuadP1>::MaxOrder() + 1; i++) {

    /* General derived parameters. */
    PetscInt num_dof_dim = i + 1;
    RealVec weights = TensorQuad<QuadP1>::GllIntegrationWeightsForOrder(i);
    RealVec points  = TensorQuad<QuadP1>::GllPointsForOrder(i);

    /* Mock constructing an element with some order. */
    PetscOptionsSetValue(NULL, "--polynomial-order", std::to_string(i).c_str());
    options->setOptions();
    TensorQuad<QuadP1> test_quad(options);
    test_quad.SetVtxCrd(vtx);
    test_quad.SetCplEdg({0, 1, 2, 3});

    /* Allocate vectors to hold approximate solutions. */
    RealMat test_field_grad(test_quad.NumIntPnt(), 2);

    /* Keep track of where we are in the tensor basis. */
    PetscInt rp = 0, sp = 0;

    /* Loop over every integration point. */
    for (PetscInt p = 0; p < test_quad.NumIntPnt(); p++) {

      /* Loop over the tensor basis. */
      for (PetscInt s = 0; s < num_dof_dim; s++) {
        for (PetscInt r = 0; r < num_dof_dim; r++) {

          /* Turn on one polynomial (dof) at a time. This is the polynomial that centered at this integration point. */
          RealVec dummy_actual = RealVec::Zero(test_quad.NumIntPnt());
          dummy_actual(r + s * num_dof_dim) = 1;

          /* Compute the gradient of this polynomial at the integration point p. */
          test_field_grad(r + s * num_dof_dim, 0) = test_quad.computeGradient(dummy_actual).col(0)(p);
          test_field_grad(r + s * num_dof_dim, 1) = test_quad.computeGradient(dummy_actual).col(1)(p);

        }
      }

      /* Get the analytical sympy generated derivative of lagrange polynomial r + s * num_dof_dim at point p. */
      RealMat analytic_grad = derivative4order(points(rp), points(sp), i);

      /* "Turn on" this the derivative belonging to integration point p. */
      RealMat test_grad_field = RealMat::Zero(test_quad.NumIntPnt(), 2);
      test_grad_field(p, 0) = 1; test_grad_field(p, 1) = 1;


      /********************************/
      /********** ASSERTIONS **********/
      /********************************/

      /* Require that the gradient of all lagrange polynomials tabulated at point p is correct. */
      REQUIRE(test_field_grad.isApprox(analytic_grad));

      /* Require that \grad test \times \grad field = 0. */
      REQUIRE(test_quad.applyGradTestAndIntegrate(test_grad_field).sum() == Approx(0.0));

      /* Require that integrals along appropriate coupling edges are correct. Since we're integrating against
       * 1, we just expect the integration weights back if an edge is hit. */
      for (int edge: {0, 1, 2, 3}) {
        PetscReal edge_val = test_quad.applyTestAndIntegrateEdge(test_grad_field.col(0), edge).sum();
        if      (sp == 0 && edge == 0) {
          REQUIRE(edge_val == Approx(weights(rp)));
        }
        else if (rp == 0 && edge == 3) {
          REQUIRE(edge_val == Approx(weights(sp)));
        }
        else if (sp == (num_dof_dim - 1) && edge == 2) {
          REQUIRE(edge_val == Approx(weights(rp)));
        }
        else if (rp == (num_dof_dim - 1) && edge == 1) {
          REQUIRE(edge_val == Approx(weights(sp)));
        }
        else {
          REQUIRE(edge_val == 0.0);
        }

      }

      /* Advance through the tensor basis. */
      rp++; if (rp == num_dof_dim) { rp = 0; sp++; }

    }

    /*
     * Interpolate of delta function. We place a delta function at the center of
     * the reference element, and ensure that, when the test functions are applied and
     * integrated against, the correct value of 1.0 is returned. This has the added
     * value of also testing applyTestAndIntegrate().
     */
    RealVec2 pnt (0.0, 0.0);
    RealVec coefficients = test_quad.getDeltaFunctionCoefficients(pnt);
    REQUIRE(test_quad.applyTestAndIntegrate(coefficients).sum() == Approx(1.0));

    /* Test that edges integrate to zero if they are set to zero (dirichlet). */
    for (int edge: {0, 1, 2, 3}) {
      RealVec test_edge = RealVec::Zero(test_quad.NumIntPnt());
      test_quad.setEdgeToValue(edge, 1.0, test_edge);
      REQUIRE(test_quad.applyTestAndIntegrateEdge(test_edge, edge).sum() == Approx(2.0));
    }

    /* TODO: Make this test tighter. */
    /* Test that the parameters interpolate. */
    RealVec4 par(1.0, 1.0, 1.0, 1.0);
    test_quad.SetVtxPar(par, "test");
    REQUIRE(test_quad.ParAtIntPts("test").sum() == Approx(test_quad.NumIntPnt()));

  }

  // Mock sources and receivers.
  PetscOptionsSetValue(NULL, "--number-of-sources", "2");
  PetscOptionsSetValue(NULL, "--source-type", "ricker");
  PetscOptionsSetValue(NULL, "--source-location-x", "50000,50000");
  PetscOptionsSetValue(NULL, "--source-location-y", "80000,90000");
  PetscOptionsSetValue(NULL, "--ricker-amplitude", "10,20");
  PetscOptionsSetValue(NULL, "--ricker-time-delay", "0.1,0.01");
  PetscOptionsSetValue(NULL, "--ricker-center-freq", "50,60");
  PetscOptionsSetValue(NULL, "--number-of-receivers", "2");
  PetscOptionsSetValue(NULL, "--receiver-names", "rec1,rec2");
  PetscOptionsSetValue(NULL, "--receiver-location-x", "50000,50000");
  PetscOptionsSetValue(NULL, "--receiver-location-y", "80000,90000");

  options->setOptions();

  // Vertex coordiantes.
  std::vector<QuadVtx> true_vtx(4);
  true_vtx[0] << 0,   0,   5e4, 0,   5e4, 8e4, 0,   8e4;
  true_vtx[1] << 5e4, 0,   1e5, 0,   1e5, 8e4, 5e4, 8e4;
  true_vtx[2] << 0,   8e4, 5e4, 8e4, 5e4, 1e5, 0,   1e5;
  true_vtx[3] << 5e4, 8e4, 1e5, 8e4, 1e5, 1e5, 5e4, 1e5;

  // True source locations.
  vector<PetscInt> locations {1, 0, 1, 0};

  auto elements = problem->initializeElements(mesh, model, options);
  std::vector<TensorQuad<QuadP1>*> p1_test;
  for (auto &e: elements) { p1_test.push_back(dynamic_cast<TensorQuad<QuadP1>*> (e.get())); }
  for (PetscInt i = 0; i < true_vtx.size(); i++) {
    REQUIRE(p1_test[i]->VtxCrd().isApprox(true_vtx[i]));
    REQUIRE(p1_test[i]->Sources().size() == locations[i]);
    REQUIRE(p1_test[i]->Receivers().size() == locations[i]);
  }

  // test for outward edge normals.
  vector<RealVec2> normals(4); vector<PetscInt> edg_elm_zero {13, 14, 15, 16};
  normals[0] << 0, -1; normals[1] << 1, 0; normals[2] << 0, 1; normals[3] << -1, 0;
  for (PetscInt i = 0; i < edg_elm_zero.size(); i++) {
    REQUIRE(p1_test[0]->getEdgeNormal(edg_elm_zero[i]).isApprox(normals[i]));
  }

}
