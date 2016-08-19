#include "catch.h"
#include <petsc.h>
#include <salvus.h>

extern "C" {
#include <Element/HyperCube/Autogen/hex_autogen.h>
}

RealMat derivative4order(const PetscReal r, const PetscReal s, const PetscReal t,
                         const PetscInt order) {

  PetscInt size = (order + 1) * (order + 1) * (order + 1);
  RealVec test_r(size), test_s(size), test_t(size);

  if (order == 1) {
    interpolate_r_derivative_order1_hex(s, t, test_r.data());
    interpolate_s_derivative_order1_hex(r, t, test_s.data());
    interpolate_t_derivative_order1_hex(r, s, test_t.data());
  } else if (order == 2) {
    interpolate_r_derivative_order2_hex(r, s, t, test_r.data());
    interpolate_s_derivative_order2_hex(r, s, t, test_s.data());
    interpolate_t_derivative_order2_hex(r, s, t, test_t.data());
  } else if (order == 3) {
    interpolate_r_derivative_order3_hex(r, s, t, test_r.data());
    interpolate_s_derivative_order3_hex(r, s, t, test_s.data());
    interpolate_t_derivative_order3_hex(r, s, t, test_t.data());
  } else if (order == 4) {
    interpolate_r_derivative_order4_hex(r, s, t, test_r.data());
    interpolate_s_derivative_order4_hex(r, s, t, test_s.data());
    interpolate_t_derivative_order4_hex(r, s, t, test_t.data());
  } else if (order == 5) {
    interpolate_r_derivative_order5_hex(r, s, t, test_r.data());
    interpolate_s_derivative_order5_hex(r, s, t, test_s.data());
    interpolate_t_derivative_order5_hex(r, s, t, test_t.data());
  } else if (order == 6) {
    interpolate_r_derivative_order6_hex(r, s, t, test_r.data());
    interpolate_s_derivative_order6_hex(r, s, t, test_s.data());
    interpolate_t_derivative_order6_hex(r, s, t, test_t.data());
  } else if (order == 7) {
    interpolate_r_derivative_order7_hex(r, s, t, test_r.data());
    interpolate_s_derivative_order7_hex(r, s, t, test_s.data());
    interpolate_t_derivative_order7_hex(r, s, t, test_t.data());
  } else if (order == 8) {
    interpolate_r_derivative_order8_hex(r, s, t, test_r.data());
    interpolate_s_derivative_order8_hex(r, s, t, test_s.data());
    interpolate_t_derivative_order8_hex(r, s, t, test_t.data());
  } else if (order == 9) {
    interpolate_r_derivative_order9_hex(r, s, t, test_r.data());
    interpolate_s_derivative_order9_hex(r, s, t, test_s.data());
    interpolate_t_derivative_order9_hex(r, s, t, test_t.data());
  } else {
    ERROR() << "Order " << order << " not supported";
  }
  RealMat ret(size, 3);
  ret.col(0) = test_r; ret.col(1) = test_s; ret.col(2) = test_t;
  return ret;
}

TEST_CASE("Test tensor hex", "[tensor_hex]") {

  /* General options. */
  PetscOptionsClear(NULL);
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--polynomial-order", "1", NULL };

  /* Fake setting via command line. */
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(NULL, &argc, &argv, NULL);

  /* Initialize dummy mesh and model. */
  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  /* Setup a reference element. */
  PetscInt num_dim = 3;
  HexVtx vtx;
  vtx <<
      -1, -1, -1,
      -1, +1, -1,
      +1, +1, -1,
      +1, -1, -1,
      -1, -1, +1,
      +1, -1, +1,
      +1, +1, +1,
      -1, +1, +1;

  /* Looping over polynomial orders < 7. Stop at 6 otherwise it takes too long! */
  for (PetscInt i = 1; i < 7; i++) {
    /* General derived parameters. */
    PetscInt num_dof_dim = i + 1;
    RealVec weights = Hexahedra<HexP1>::GllIntegrationWeights(i);
    RealVec points  = Hexahedra<HexP1>::GllPointsForOrder(i);

    /* Construct an element with some polynomial order. */
    PetscOptionsSetValue(NULL, "--polynomial-order", std::to_string(i).c_str());
    options->setOptions();
    Hexahedra<HexP1> test_hex(options);
    test_hex.SetVtxCrd(vtx);

    /* Allocate vectors to hold approximate solutions. */
    RealMat test_field_grad(test_hex.NumIntPnt(), num_dim);

    /* Keep track of where we are in the tensor basis. */
    PetscInt rp = 0, sp = 0, tp = 0;

    /* Loop over every integration point. */
    for (PetscInt p = 0; p < test_hex.NumIntPnt(); p++) {

      /* Loop over tensor basis. */
      for (PetscInt t = 0; t < num_dof_dim; t++) {
        for (PetscInt s = 0; s < num_dof_dim; s++) {
          for (PetscInt r = 0; r < num_dof_dim; r++) {

            PetscInt ind = r + s * num_dof_dim + t * num_dof_dim * num_dof_dim;
            /* Turn on one polynomial (dof) at a time.
             * This is the polynomial that centered at this integration point. */
            RealVec dummy_actual = RealVec::Zero(test_hex.NumIntPnt());
            dummy_actual(ind) = 1;

            /* Compute the gradient of this polynomial at the integration point p. */
            test_field_grad(ind, 0) = test_hex.computeGradient(dummy_actual).col(0)(p);
            test_field_grad(ind, 1) = test_hex.computeGradient(dummy_actual).col(1)(p);
            test_field_grad(ind, 2) = test_hex.computeGradient(dummy_actual).col(2)(p);


          }
        }
      }

      /* Get the analytical sympy generated derivative of lagrange polynomial r + s * num_dof_dim at point p. */
      RealMat analytic_grad = derivative4order(points(rp), points(sp), points(tp), i);

      /* "Turn on" this the derivative belonging to integration point p. */
      RealMat test_grad_field = RealMat::Zero(test_hex.NumIntPnt(), 3);
      test_grad_field(p, 0) = 1; test_grad_field(p, 1) = 1; test_grad_field(p, 2) = 1;

      /* Require that \grad test \times \grad field = 0. */
      REQUIRE(test_hex.applyGradTestAndIntegrate(test_grad_field).sum() == Approx(0.0));

      /********************************/
      /********** ASSERTIONS **********/
      /********************************/

      /* Require that the gradient of all lagrange polynomials tabulated at point p is correct. */
      REQUIRE(test_field_grad.isApprox(analytic_grad));

      for (int edge: {0, 1, 2, 3, 4, 5}) {
        PetscReal edge_val = test_hex.applyTestAndIntegrateEdge(test_grad_field.col(0), edge).sum();
        switch (edge) {
          case 0: /* On bottom face. */
            if (tp == 0) {
              REQUIRE(edge_val == Approx(weights(rp) * weights(sp)));
            } else {
              REQUIRE(edge_val == Approx(0.0));
            }
            break;
          case 1: /* On top face. */
            if (tp == num_dof_dim - 1) {
              REQUIRE(edge_val == Approx(weights(rp) * weights(sp)));
            } else {
              REQUIRE(edge_val == Approx(0.0));
            }
            break;
          case 2: /* On left face. */
            if (sp == 0) {
              REQUIRE(edge_val == Approx(weights(rp) * weights(tp)));
            } else {
              REQUIRE(edge_val == Approx(0.0));
            }
            break;
          case 3: /* On right face. */
            if (sp == num_dof_dim - 1) {
              REQUIRE(edge_val == Approx(weights(rp) * weights(tp)));
            } else {
              REQUIRE(edge_val == Approx(0.0));
            }
            break;
          case 4: /* On back face. */
            if (rp == num_dof_dim - 1) {
              REQUIRE(edge_val == Approx(weights(sp) * weights(tp)));
            } else {
              REQUIRE(edge_val == Approx(0.0));
            }
            break;
          case 5: /* On front face. */
            if (rp == 0) {
              REQUIRE(edge_val == Approx(weights(sp) * weights(tp)));
            } else {
              REQUIRE(edge_val == Approx(0.0));
            }
            break;
        }
      }

      /* TODO: CHECK THIS. Are r and s flipped? Need to look into this. */
      /* Advance through the tensor basis. */
      rp++;
      if (rp == num_dof_dim) { rp = 0; sp++; }
      if (sp == num_dof_dim) { sp = 0; tp++; }

    }

    /*
     * Interpolate of delta function. We place a delta function at the center of
     * the reference element, and ensure that, when the test functions are applied and
     * integrated against, the correct value of 1.0 is returned. This has the added
     * value of also testing applyTestAndIntegrate().
     */
    RealVec3 pnt (0.0, 0.0, 0.0);
    RealVec coefficients = test_hex.getDeltaFunctionCoefficients(pnt);
    REQUIRE(test_hex.applyTestAndIntegrate(coefficients).sum() == Approx(1.0));

    /* Test that faces are set properly. */
    for (int edge: {0, 1, 2, 3, 4, 5}) {
      RealVec test_face = RealVec::Zero(test_hex.NumIntPnt());
      test_hex.setEdgeToValue(edge, 1.0, test_face);
      REQUIRE(test_hex.applyTestAndIntegrateEdge(test_face, edge).sum() == Approx(4.0));
    }

    /* TODO: MAKE THIS TEST TIGHTER. */
    /* Test that the parameters interpolate. */
    RealVec par(8);
    par.setConstant(1.0);
    test_hex.SetVtxPar(par, "test");
    REQUIRE(test_hex.ParAtIntPts("test").sum() == Approx(test_hex.NumIntPnt()));


  }

  // Mock sources and receivers.
  std::string e_file = "small_hex_mesh_to_test_sources.e";
  PetscOptionsSetValue(NULL, "--mesh-file", e_file.c_str());
  PetscOptionsSetValue(NULL, "--model-file", e_file.c_str());
  PetscOptionsSetValue(NULL, "--number-of-sources", "2");
  PetscOptionsSetValue(NULL, "--source-type", "ricker");
  PetscOptionsSetValue(NULL, "--source-location-x", "50000,50000");
  PetscOptionsSetValue(NULL, "--source-location-y", "50000,90000");
  PetscOptionsSetValue(NULL, "--source-location-z", "50000,90000");
  PetscOptionsSetValue(NULL, "--source-num-components", "3,3");
  PetscOptionsSetValue(NULL, "--ricker-amplitude", "1,1,1");
  PetscOptionsSetValue(NULL, "--ricker-time-delay", "1,1,1");
  PetscOptionsSetValue(NULL, "--ricker-center-freq", "1,1,1");
  PetscOptionsSetValue(NULL, "--receiver-file-name", "mock.h5");
  PetscOptionsSetValue(NULL, "--number-of-receivers", "2");
  PetscOptionsSetValue(NULL, "--receiver-names", "rec1,rec2");
  PetscOptionsSetValue(NULL, "--receiver-location-x", "50000,50000");
  PetscOptionsSetValue(NULL, "--receiver-location-y", "50000,90000");
  PetscOptionsSetValue(NULL, "--receiver-location-z", "50000,90000");

  /* TODO: Fix this. */
  options->SetDimension(3);
  options->setOptions();

  /* Set up a simple mesh. */
  std::unique_ptr<Problem>      problem(Problem::Factory(options));
  std::unique_ptr<ExodusModel>  model(new ExodusModel(options));
  std::unique_ptr<Mesh>         mesh(new Mesh(options));

  mesh->read();
  model->read();
  mesh->setupTopology(model, options);
  auto elements = problem->initializeElements(mesh, model, options);
  mesh->setupGlobalDof(elements[0], options);

  /* Known vertex coordinates. */
  std::vector<PetscInt> locations {1, 0, 0, 0, 0, 0, 1, 0};
  std::vector<HexVtx> true_vtx(8);
  true_vtx[0] <<
    0, 0, 0, 0, 0, 50000, 0, 50000, 50000, 0, 50000, 0, 50000, 0, 0, 50000,
    50000, 0, 50000, 50000, 50000, 50000, 0, 50000;
  true_vtx[1] << 50000, 0, 0, 50000, 0, 50000, 50000, 50000, 50000, 50000,
      50000, 0, 100000, 0, 0, 100000, 50000, 0, 100000, 50000, 50000,
      100000, 0, 50000;
  true_vtx[2] << 0, 50000, 0, 0, 50000, 50000, 0, 100000, 50000, 0, 100000,
      0, 50000, 50000, 0, 50000, 100000, 0, 50000, 100000, 50000, 50000,
      50000, 50000;
  true_vtx[3] << 50000, 50000, 0, 50000, 50000, 50000, 50000, 100000,
      50000, 50000, 100000, 0, 100000, 50000, 0, 100000, 100000, 0,
      100000, 100000, 50000, 100000, 50000, 50000;
  true_vtx[4] << 0, 0, 50000, 0, 0, 100000, 0, 50000, 100000, 0, 50000, 50000,
      50000, 0, 50000, 50000, 50000, 50000, 50000, 50000, 100000, 50000, 0,
      100000;
  true_vtx[5] << 50000, 0, 50000, 50000, 0, 100000, 50000, 50000, 100000,
      50000, 50000, 50000, 100000, 0, 50000, 100000, 50000, 50000, 100000,
      50000, 100000, 100000, 0, 100000;
  true_vtx[6] << 0, 50000, 50000, 0, 50000, 100000, 0, 100000, 100000, 0,
      100000, 50000, 50000, 50000, 50000, 50000, 100000, 50000, 50000,
      100000, 100000, 50000, 50000, 100000;
  true_vtx[7] << 50000, 50000, 50000, 50000, 50000, 100000, 50000, 100000,
      100000, 50000, 100000, 50000, 100000, 50000, 50000, 100000, 100000,
      50000, 100000, 100000, 100000, 100000, 50000, 100000;

  /* Get raw elements and cast to test. */
  std::vector<Hexahedra<HexP1>*> p1_test;
  for (auto &e: elements) { p1_test.push_back(dynamic_cast<Hexahedra<HexP1>*> (e.get())); }
  for (PetscInt i = 0; i < true_vtx.size(); i++) {
    REQUIRE(p1_test[i]->VtxCrd().isApprox(true_vtx[i]));
    REQUIRE(p1_test[i]->Sources().size() == locations[i]);
    REQUIRE(p1_test[i]->Receivers().size() == locations[i]);
  }

  /* Require that proper errors are thrown. */
  PetscOptionsSetValue(NULL, "--polynomial-order", "10");
  options->setOptions();
  REQUIRE_THROWS_AS(new Hexahedra<HexP1>(options), std::runtime_error);
  REQUIRE_THROWS_AS(Hexahedra<HexP1>::GllPointsForOrder(
      Hexahedra<HexP1>::MaxOrder()+1), std::runtime_error);
  REQUIRE_THROWS_AS(Hexahedra<HexP1>::GllIntegrationWeights(
      Hexahedra<HexP1>::MaxOrder()+1), std::runtime_error);
  REQUIRE_THROWS_AS(Hexahedra<HexP1>::setupGradientOperator(
      Hexahedra<HexP1>::MaxOrder()+1), std::runtime_error);
  REQUIRE_THROWS_AS(Hexahedra<HexP1>::interpolateLagrangePolynomials(
      0, 0, 0, Hexahedra<HexP1>::MaxOrder()+1), std::runtime_error);

}
