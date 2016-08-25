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
  }
  #if HEX_MAX_ORDER > 7
  else if (order == 8) {
    interpolate_r_derivative_order8_hex(r, s, t, test_r.data());
    interpolate_s_derivative_order8_hex(r, s, t, test_s.data());
    interpolate_t_derivative_order8_hex(r, s, t, test_t.data());
  } else if (order == 9) {
    interpolate_r_derivative_order9_hex(r, s, t, test_r.data());
    interpolate_s_derivative_order9_hex(r, s, t, test_s.data());
    interpolate_t_derivative_order9_hex(r, s, t, test_t.data());
  }
  #endif
  else {
    ERROR() << "Order " << order << " not supported";
  }
  RealMat ret(size, 3);
  ret.col(0) = test_r; ret.col(1) = test_s; ret.col(2) = test_t;
  return ret;
}

TEST_CASE("Test tensor hex", "[tensor_hex]") {

  SECTION("Mathematical operations") {

    /* General options. */
    PetscOptionsClear(NULL);
    const char *arg[] = {"salvus_test", "--testing", "true", "--polynomial-order", "1", NULL};

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
    vtx
        << -1, -1, -1, -1, +1, -1, +1, +1, -1, +1, -1, -1, -1, -1, +1, +1, -1, +1, +1, +1, +1, -1, +1, +1;

    /* Looping over polynomial orders < 7. Stop at 6 otherwise it takes too long! */
#ifdef NDEBUG    
    for (PetscInt i = 1; i < 7; i++) {
#else
      for (PetscInt i = 3; i < 3 + 1; i++) {
#endif      
      /* General derived parameters. */
      PetscInt num_dof_dim = i + 1;
      RealVec weights = Hexahedra<HexP1>::GllIntegrationWeights(i);
      RealVec points = Hexahedra<HexP1>::GllPointsForOrder(i);

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
        test_grad_field(p, 0) = 1;
        test_grad_field(p, 1) = 1;
        test_grad_field(p, 2) = 1;

        /* Require that \grad test \times \grad field = 0. */
        REQUIRE(test_hex.applyGradTestAndIntegrate(test_grad_field).sum() == Approx(0.0));

        /********************************/
        /********** ASSERTIONS **********/
        /********************************/

        /* Require that the gradient of all lagrange polynomials tabulated at point p is correct. */
        REQUIRE(test_field_grad.isApprox(analytic_grad));

        for (int edge: {0, 1, 2, 3, 4, 5}) {
          PetscReal
              edge_val = test_hex.applyTestAndIntegrateEdge(test_grad_field.col(0), edge).sum();
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
        if (rp == num_dof_dim) {
          rp = 0;
          sp++;
        }
        if (sp == num_dof_dim) {
          sp = 0;
          tp++;
        }

      }

      /*
       * Interpolate of delta function. We place a delta function at the center of
       * the reference element, and ensure that, when the test functions are applied and
       * integrated against, the correct value of 1.0 is returned. This has the added
       * value of also testing applyTestAndIntegrate().
       */
      RealVec3 pnt(0.0, 0.0, 0.0);
      RealVec coefficients = test_hex.getDeltaFunctionCoefficients(pnt);
      REQUIRE(test_hex.applyTestAndIntegrate(coefficients).sum() == Approx(1.0));

      /* Test that faces are set properly. */
      for (int edge: {0, 1, 2, 3, 4, 5}) {
        RealVec test_face = RealVec::Zero(test_hex.NumIntPnt());
        for (PetscInt j: test_hex.getDofsOnFace(edge)) { test_face(j) = 1.0; }
//      test_hex.setFaceToValue(edge, 1.0, test_face);
        REQUIRE(test_hex.applyTestAndIntegrateEdge(test_face, edge).sum() == Approx(4.0));
      }

      /* TODO: MAKE THIS TEST TIGHTER. */
      /* Test that the parameters interpolate. */
      RealVec par(8);
      par.setConstant(1.0);
      test_hex.SetVtxPar(par, "test");
      REQUIRE(test_hex.ParAtIntPts("test").sum() == Approx(test_hex.NumIntPnt()));

    }
  }

  SECTION("Test closure for order 4") {

    PetscOptionsSetValue(NULL, "--polynomial-order", "4");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();
    Hexahedra<HexP1> test_hex(options);
    SECTION("Vertices") {
      REQUIRE(test_hex.getDofsOnVtx(0) == 0);
      REQUIRE(test_hex.getDofsOnVtx(1) == 20);
      REQUIRE(test_hex.getDofsOnVtx(2) == 24);
      REQUIRE(test_hex.getDofsOnVtx(3) == 4);
      REQUIRE(test_hex.getDofsOnVtx(4) == 100);
      REQUIRE(test_hex.getDofsOnVtx(5) == 104);
      REQUIRE(test_hex.getDofsOnVtx(6) == 124);
      REQUIRE(test_hex.getDofsOnVtx(7) == 120);
      REQUIRE_THROWS_AS(test_hex.getDofsOnVtx(8), std::runtime_error);
    }
    SECTION("Edges") {
      std::vector<PetscInt> e0  = {0,     5,  10,  15,  20};
      std::vector<PetscInt> e1  = {20,   21,  22,  23,  24};
      std::vector<PetscInt> e2  = {4,     9,  14,  19,  24};
      std::vector<PetscInt> e3  = {0,     1,   2,   3,   4};
      std::vector<PetscInt> e4  = {100, 101, 102, 103, 104};
      std::vector<PetscInt> e5  = {104, 109, 114, 119, 124};
      std::vector<PetscInt> e6  = {120, 121, 122, 123, 124};
      std::vector<PetscInt> e7  = {100, 105, 110, 115, 120};
      std::vector<PetscInt> e8  = {4,    29,  54,  79, 104};
      std::vector<PetscInt> e9  = {0,    25,  50,  75, 100};
      std::vector<PetscInt> e10 = {20,   45,  70,  95, 120};
      std::vector<PetscInt> e11 = {24,   49,  74,  99, 124};
      REQUIRE(test_hex.getDofsOnEdge(0)  == e0);
      REQUIRE(test_hex.getDofsOnEdge(1)  == e1);
      REQUIRE(test_hex.getDofsOnEdge(2)  == e2);
      REQUIRE(test_hex.getDofsOnEdge(3)  == e3);
      REQUIRE(test_hex.getDofsOnEdge(4)  == e4);
      REQUIRE(test_hex.getDofsOnEdge(5)  == e5);
      REQUIRE(test_hex.getDofsOnEdge(6)  == e6);
      REQUIRE(test_hex.getDofsOnEdge(7)  == e7);
      REQUIRE(test_hex.getDofsOnEdge(8)  == e8);
      REQUIRE(test_hex.getDofsOnEdge(9)  == e9);
      REQUIRE(test_hex.getDofsOnEdge(10) == e10);
      REQUIRE(test_hex.getDofsOnEdge(11) == e11);
      REQUIRE_THROWS_AS(test_hex.getDofsOnEdge(12), std::runtime_error);
    }

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

  SECTION("Integration with a simple mesh") {

    // Mock sources and receivers.
    std::string e_file = "small_hex_mesh_to_test_sources.e";
    PetscOptionsSetValue(NULL, "--mesh-file", e_file.c_str());
    PetscOptionsSetValue(NULL, "--model-file", e_file.c_str());
    PetscOptionsSetValue(NULL, "--number-of-sources", "2");
    PetscOptionsSetValue(NULL, "--source-type", "ricker");
    PetscOptionsSetValue(NULL, "--source-location-x", "50000,50000");
    PetscOptionsSetValue(NULL, "--source-location-y", "50000,90000");
    PetscOptionsSetValue(NULL, "--source-location-z", "50000,90000");
    PetscOptionsSetValue(NULL, "--ricker-amplitude", "1,1,1");
    PetscOptionsSetValue(NULL, "--ricker-time-delay", "1,1,1");
    PetscOptionsSetValue(NULL, "--ricker-center-freq", "1,1,1");
    PetscOptionsSetValue(NULL, "--receiver-file-name", "mock.h5");
    PetscOptionsSetValue(NULL, "--number-of-receivers", "2");
    PetscOptionsSetValue(NULL, "--receiver-names", "rec1,rec2");
    PetscOptionsSetValue(NULL, "--receiver-location-x", "50000,50000");
    PetscOptionsSetValue(NULL, "--receiver-location-y", "50000,90000");
    PetscOptionsSetValue(NULL, "--receiver-location-z", "50000,90000");
    PetscOptionsSetValue(NULL, "--polynomial-order", "5");

    /* TODO: Fix this. */
    std::unique_ptr<Options> options(new Options);
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
}

TEST_CASE("test closure mapping","[hex/closure]") {

  std::string e_file = "hex_eigenfunction.e";

  PetscOptionsClear(NULL);
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--mesh-file", e_file.c_str(),
      "--model-file", e_file.c_str(),
      "--time-step", "1e-2",
      "--polynomial-order", "3", NULL};
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(NULL, &argc, &argv, NULL);

  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  std::unique_ptr<Problem> problem(Problem::Factory(options));
  std::unique_ptr<ExodusModel> model(new ExodusModel(options));
  std::unique_ptr<Mesh> mesh(Mesh::Factory(options));

  model->read();
  PetscInt cells[] =
      {0, 1, 2, 3,
       4, 5, 6, 7,

       12, 13, 14, 15, // element below
       0, 3, 2, 1,

       4, 7, 6, 5,    // element above
       8, 9, 10, 11,

       16, 0, 3, 19,  // element front
       17, 18, 5, 4,

       1, 20, 23, 2,  // element back
       7, 6, 22, 21,

       3, 2, 25, 24,  // element right
       5, 27, 26, 6,

       28, 29, 1, 0,  // element left
       31, 4, 7, 30 };

  // 4,7,6,5, // bottom o:-4
  // 8,9,10,11 // top o:0
  // 5,4,7,6, // bottom o:-1
  // 9,10,11,8 // top o:0
  // 6,5,4,7, // bottom o:-2
  // 10,11,8,9 // top o:0
  // 7,6,5,4, // bottom o:-3
  // 10,8,9,11 // top o:0

  PetscReal vertex_coords[] =
      {+0.0, +0.0, +0.0, //0
       +0.0, +1.0, +0.0,
       +1.0, +1.0, +0.0,
       +1.0, +0.0, +0.0,

       +0.0, +0.0, +1.0, // 4
       +1.0, +0.0, +1.0,
       +1.0, +1.0, +1.0,
       +0.0, +1.0, +1.0,

       +0.0, +0.0, +2.0, // 8
       +1.0, +0.0, +2.0,
       +1.0, +1.0, +2.0,
       +0.0, +1.0, +2.0,

       +0.0, +0.0, -1.0, // 12
       +0.0, +1.0, -1.0,
       +1.0, +1.0, -1.0,
       +1.0, +0.0, -1.0,

       +0.0, -1.0, +0.0, // 16
       +0.0, -1.0, +1.0,
       +1.0, -1.0, +1.0,
       +1.0, -1.0, +0.0,

       +0.0, +2.0, +0.0, // 20
       +0.0, +2.0, +1.0,
       +1.0, +2.0, +1.0,
       +1.0, +2.0, +0.0,

       +2.0, +0.0, +0.0, // 24
       +2.0, +1.0, +0.0,
       +2.0, +1.0, +1.0,
       +2.0, +0.0, +1.0,

       -1.0, +0.0, +0.0, // 28
       -1.0, +1.0, +0.0,
       -1.0, +1.0, +1.0,
       -1.0, +0.0, +1.0 };

  /* Creat a mesh with the vertices defined above. */
  PetscInt dimen = 3, num_cells = 7, num_verts = 32, verts_per_cell = 8;
  mesh->read(dimen, num_cells, num_verts, verts_per_cell, cells, vertex_coords);

  /* Setup topology from model and mesh. */
  mesh->setupTopology(model, options);

  /* Setup elements from model and topology. */
  auto elements = problem->initializeElements(mesh, model, options);

  /* Setup global degrees of freedom based on element 0. */
  mesh->setupGlobalDof(elements[0], options);

  auto fields = problem->initializeGlobalDofs(elements, mesh);
  auto test_hex = dynamic_cast<Hexahedra<HexP1>*> (elements[0].get());

  SECTION("Ensure face closure is set properly") {

    std::vector<PetscInt> insert_element {1, 2, 3, 4, 5, 6};
    std::vector<PetscInt> insert_faces   {1, 0, 3, 2, 5, 4};

    RealVec un = RealVec::Zero(elements[0]->NumIntPnt());

    RealVec bot(un.size()), top(un.size()), lft(un.size()), rgt(un.size()),
        fnt(un.size()), bck(un.size());
    bot << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0;
    top << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1;
    fnt << 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0;
    bck << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1;
    rgt << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
        1, 0, 0, 0, 1;
    lft << 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0,
        0, 1, 0, 0, 0;

    for (PetscInt i = 0; i < insert_element.size(); i++) {

      /* Insert a vector of ones of some face of a neighbouring element. */
      for (PetscInt j: test_hex->getDofsOnFace(insert_faces[i])) { un(j) = 1.0; }
      problem->insertElementalFieldIntoMesh("u", insert_element[i], elements[0]->ClsMap(), un,
                                            mesh->DistributedMesh(), mesh->MeshSection(),
                                            fields);

      switch (insert_element[i]) {

        /* Get field on central element (number 0), one face of which should be set
         * to one (others to zero). */
        case 1:
          REQUIRE(problem->getFieldOnElement("u", 0, elements[0]->ClsMap(),
                                             mesh->DistributedMesh(),
                                             mesh->MeshSection(), fields).isApprox(bot));
          break;
        case 2:
          REQUIRE(problem->getFieldOnElement("u", 0, elements[0]->ClsMap(),
                                             mesh->DistributedMesh(),
                                             mesh->MeshSection(), fields).isApprox(top));
          break;
        case 3:
          REQUIRE(problem->getFieldOnElement("u", 0, elements[0]->ClsMap(),
                                             mesh->DistributedMesh(),
                                             mesh->MeshSection(), fields).isApprox(fnt));
          break;
        case 4:
          REQUIRE(problem->getFieldOnElement("u", 0, elements[0]->ClsMap(),
                                             mesh->DistributedMesh(),
                                             mesh->MeshSection(), fields).isApprox(bck));
          break;
        case 5:
          REQUIRE(problem->getFieldOnElement("u", 0, elements[0]->ClsMap(),
                                             mesh->DistributedMesh(),
                                             mesh->MeshSection(), fields).isApprox(rgt));
          break;
        case 6:
          REQUIRE(problem->getFieldOnElement("u", 0, elements[0]->ClsMap(),
                                             mesh->DistributedMesh(),
                                             mesh->MeshSection(), fields).isApprox(lft));
        default:
          break;

      }

      /* Reset that face to zero. */
      for (PetscInt j: test_hex->getDofsOnFace(insert_faces[i])) { un(j) = 0.0; }
      problem->insertElementalFieldIntoMesh("u", insert_element[i], elements[0]->ClsMap(), un,
                                            mesh->DistributedMesh(), mesh->MeshSection(),
                                            fields);
    }
  }

  SECTION("Closure routines throw on bad values") {

    RealVec dum(1);
    REQUIRE_THROWS_AS(test_hex->applyTestAndIntegrateEdge(dum, 6), std::runtime_error);
    REQUIRE_THROWS_AS(test_hex->getDofsOnFace(6), std::runtime_error);

  }


}

