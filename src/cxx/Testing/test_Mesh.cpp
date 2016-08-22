#include <iostream>
#include <salvus.h>
#include "catch.h"


TEST_CASE("Unit test mesh", "[mesh]") {

  PetscOptionsClear(NULL);
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--mesh-file", "quad_eigenfunction",
      "--polynomial-order", "1", NULL};

  /* Fake setting via command line. */
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(NULL, &argc, &argv, NULL);

  SECTION("Fail with non-existant mesh file.") {
    /* Start with the wrong mesh name. */
    std::unique_ptr<Options> options(new Options());
    options->setOptions();
    auto mesh = Mesh::Factory(options);
    REQUIRE_THROWS_AS(mesh->read(), std::runtime_error);
  }

  SECTION("Fail when we forget to call read()") {

    std::unique_ptr<Options> options(new Options());
    options->setOptions();
    std::unique_ptr<ExodusModel> model(new ExodusModel(options));
    auto mesh = Mesh::Factory(options);

    REQUIRE_THROWS_AS(mesh->setupTopology(model, options), std::runtime_error);

  }

  SECTION("Fail when we try to setup global dofs before we set up topology") {

    PetscOptionsSetValue(NULL, "--mesh-file", "quad_eigenfunction.e");
    PetscOptionsSetValue(NULL, "--model-file", "quad_eigenfunction.e");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    auto mesh = Mesh::Factory(options);
    mesh->read();

    REQUIRE_THROWS_AS(mesh->setupGlobalDof(Element::Factory("quad", {"fluid"}, {}, options),
                                           options), std::runtime_error);


  }

  SECTION("Correctly initialize DM (2D quad).") {

    PetscOptionsSetValue(NULL, "--mesh-file", "quad_eigenfunction.e");
    PetscOptionsSetValue(NULL, "--model-file", "quad_eigenfunction.e");
    PetscOptionsSetValue(NULL, "--homogeneous-dirichlet", "x0");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    /* Ensure mesh was read. */
    auto mesh = Mesh::Factory(options);
    mesh->read();
    REQUIRE(mesh->DistributedMesh() != NULL);

    /* Attach model and global dofs. */
    std::unique_ptr<ExodusModel> model(new ExodusModel(options));
    model->read();
    mesh->setupTopology(model, options);

    /* Generate dofs from elements. */
    std::unique_ptr<Problem> problem(Problem::Factory(options));
    auto elements = problem->initializeElements(mesh, model, options);
    mesh->setupGlobalDof(elements[0], options);
    REQUIRE(mesh->MeshSection() != NULL);

    /* Ensure vertices are correct. */
    QuadVtx vtx; vtx << 0, 0, 25000, 0, 25000, 25000, 0, 25000;
    REQUIRE(mesh->getElementCoordinateClosure(0).isApprox(vtx));

    /* Ensure edges are consistent. */
    std::vector<PetscInt> true_edges = { 41, 42, 43, 44 };
    REQUIRE(mesh->EdgeNumbers(0) == true_edges);

    /* Ensure neighbouring element are consistent. */
    REQUIRE_THROWS_AS(mesh->GetNeighbouringElement(41, 0), std::runtime_error);
    REQUIRE_THROWS_AS(mesh->GetNeighbouringElement(44, 0), std::runtime_error);
    REQUIRE(mesh->GetNeighbouringElement(42, 0) == 1);
    REQUIRE(mesh->GetNeighbouringElement(43, 0) == 4);

    /* Ensure that we only return the physics of our neighbours. */
    REQUIRE(mesh->CouplingFields((0)).size() == 2);
    REQUIRE(mesh->CouplingFields((1)).size() == 3);
    REQUIRE(mesh->CouplingFields((5)).size() == 4);

    /* Ensure that we can get the proper physics from neighbour elements. */
    for (auto &f: mesh->CouplingFields(0)) {
      REQUIRE(std::get<1>(f).size() == 1);
      REQUIRE(std::get<1>(f)[0] == "fluid");
    }

    /* Ensure that we get the proper physics for our element. */
    REQUIRE(mesh->ElementFields(0)[0] == "fluid");

    /* Ensure that, if we're on a boundary, we can detect that. */
    REQUIRE(mesh->TotalCouplingFields(0)[0] == "boundary_homo_dirichlet");
    /* Only a boundary on the left edge, so element 1 should have free surface (no coupling). */
    REQUIRE(mesh->TotalCouplingFields(1).empty());

    /* Ensure that if we're in the middle of a homogeneous section, there's no coupling. */
    REQUIRE(mesh->TotalCouplingFields(5).empty());

    /* Make sure we know we're a quad. */
    REQUIRE(mesh->baseElementType() == "quad");

    /* Make sure we have the correct number of elements. */
    REQUIRE(mesh->NumberElementsLocal() == 16);

    /* Ensure that we get a correct listing of all the boundary points. */
    std::set<std::tuple<PetscInt,PetscInt>> all_boundaries_true {
        std::make_tuple(0, 16), std::make_tuple(0, 17), std::make_tuple(0, 18),
        std::make_tuple(0, 19), std::make_tuple(0, 20), std::make_tuple(0, 44),
        std::make_tuple(0, 56), std::make_tuple(0, 65), std::make_tuple(0, 74),
        std::make_tuple(1, 36), std::make_tuple(1, 37), std::make_tuple(1, 38),
        std::make_tuple(1, 39), std::make_tuple(1, 40), std::make_tuple(1, 52),
        std::make_tuple(1, 61), std::make_tuple(1, 70), std::make_tuple(1, 79),
        std::make_tuple(2, 16), std::make_tuple(2, 21), std::make_tuple(2, 26),
        std::make_tuple(2, 31), std::make_tuple(2, 36), std::make_tuple(2, 41),
        std::make_tuple(2, 45), std::make_tuple(2, 48), std::make_tuple(2, 51),
        std::make_tuple(3, 20), std::make_tuple(3, 25), std::make_tuple(3, 30),
        std::make_tuple(3, 35), std::make_tuple(3, 40), std::make_tuple(3, 73),
        std::make_tuple(3, 76), std::make_tuple(3, 78), std::make_tuple(3, 80) };

    REQUIRE(mesh->BoundaryPoints() == all_boundaries_true);

  }

  SECTION("Correctly initialize DM (3D hex).") {
    PetscOptionsSetValue(NULL, "--mesh-file",  "small_hex_mesh_to_test_sources.e");
    PetscOptionsSetValue(NULL, "--model-file", "small_hex_mesh_to_test_sources.e");
    PetscOptionsSetValue(NULL, "--homogeneous-dirichlet", "x0");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    /* Ensure mesh was read. */
    auto mesh = Mesh::Factory(options);
    mesh->read();
    REQUIRE(mesh->DistributedMesh() != NULL);

    /* Attach model and global dofs. */
    std::unique_ptr<ExodusModel> model(new ExodusModel(options));
    model->read();
    mesh->setupTopology(model, options);
    std::unique_ptr<Problem> problem(Problem::Factory(options));
    auto elements = problem->initializeElements(mesh, model, options);
    mesh->setupGlobalDof(elements[0], options);
    REQUIRE(mesh->MeshSection() != NULL);

    /* Ensure vertices are correct. */
    HexVtx vtx;
    vtx <<
        0,     0,     0,
        0,     0,     50000,
        0,     50000, 50000,
        0,     50000, 0,
        50000, 0,     0,
        50000, 50000, 0,
        50000, 50000, 50000,
        50000, 0,     50000;
    REQUIRE(mesh->getElementCoordinateClosure(0).isApprox(vtx));

    /* Ensure edges are consistent. */
    std::vector<PetscInt> true_edges = { 35, 36, 37, 38, 39, 40 };
    REQUIRE(mesh->EdgeNumbers(0) == true_edges);

    /* Make sure our neighbours are properly known. */
    REQUIRE_THROWS_AS(mesh->GetNeighbouringElement(35, 0), std::runtime_error);
    REQUIRE_THROWS_AS(mesh->GetNeighbouringElement(37, 0), std::runtime_error);
    REQUIRE_THROWS_AS(mesh->GetNeighbouringElement(40, 0), std::runtime_error);
    REQUIRE(mesh->GetNeighbouringElement(36, 0) == 1);
    REQUIRE(mesh->GetNeighbouringElement(38, 0) == 4);
    REQUIRE(mesh->GetNeighbouringElement(39, 0) == 2);

    /* TODO: Modify the mesh to include elements not on the edges. */
    /* Ensure that we only return the physics of our neighbours. */
    REQUIRE(mesh->CouplingFields((0)).size() == 3);
    REQUIRE(mesh->CouplingFields((1)).size() == 3);
    REQUIRE(mesh->CouplingFields((5)).size() == 3);

    /* Ensure that we can get the proper physics from neighbour elements. */
    for (auto &f: mesh->CouplingFields(0)) {
      REQUIRE(std::get<1>(f).size() == 1);
      REQUIRE(std::get<1>(f)[0] == "fluid");
    }

    /* Ensure that we get the proper physics for our element. */
    REQUIRE(mesh->ElementFields(0)[0] == "fluid");

    /* Ensure that, if we're on a boundary, we can detect that. */
    REQUIRE(mesh->TotalCouplingFields(0)[0] == "boundary_homo_dirichlet");
    /* We've only set x0 to a boundary, so all other interfaces should be empty. */
    REQUIRE(mesh->TotalCouplingFields(1).empty());

    /* TODO: See above TODO. */
    /* Ensure that if we're in the middle of a homogeneous section, there's no coupling. */
    /* REQUIRE(mesh->TotalCouplingFields(5).empty() == true); */

    /* Ensure we know we're a hex. */
    REQUIRE(mesh->baseElementType() == "hex");

    /* Ensure we have the proper number of elements. */
    REQUIRE(mesh->NumberElementsLocal() == 8);

    /* Ensure that we get a correct listing of all the boundary points. */

    std::set<std::tuple<PetscInt,PetscInt>> all_boundaries_true {
        std::make_tuple(0, 8), std::make_tuple(0, 9), std::make_tuple(0, 10),
        std::make_tuple(0, 11), std::make_tuple(0, 12), std::make_tuple(0, 13),
        std::make_tuple(0, 14), std::make_tuple(0, 15), std::make_tuple(0, 16),
        std::make_tuple(0, 35), std::make_tuple(0, 46), std::make_tuple(0, 55),
        std::make_tuple(0, 64), std::make_tuple(0, 71), std::make_tuple(0, 72),
        std::make_tuple(0, 73), std::make_tuple(0, 74), std::make_tuple(0, 91),
        std::make_tuple(0, 92), std::make_tuple(0, 93), std::make_tuple(0, 104),
        std::make_tuple(0, 105), std::make_tuple(0, 106), std::make_tuple(0, 117),
        std::make_tuple(0, 118), std::make_tuple(1, 26), std::make_tuple(1, 27),
        std::make_tuple(1, 28), std::make_tuple(1, 29), std::make_tuple(1, 30),
        std::make_tuple(1, 31), std::make_tuple(1, 32), std::make_tuple(1, 33),
        std::make_tuple(1, 34), std::make_tuple(1, 41), std::make_tuple(1, 51),
        std::make_tuple(1, 60), std::make_tuple(1, 68), std::make_tuple(1, 83),
        std::make_tuple(1, 84), std::make_tuple(1, 85), std::make_tuple(1, 86),
        std::make_tuple(1, 99), std::make_tuple(1, 100), std::make_tuple(1, 101),
        std::make_tuple(1, 112), std::make_tuple(1, 113), std::make_tuple(1, 114),
        std::make_tuple(1, 122), std::make_tuple(1, 123), std::make_tuple(2, 8),
        std::make_tuple(2, 9), std::make_tuple(2, 10), std::make_tuple(2, 17),
        std::make_tuple(2, 18), std::make_tuple(2, 19), std::make_tuple(2, 26),
        std::make_tuple(2, 27), std::make_tuple(2, 28), std::make_tuple(2, 40),
        std::make_tuple(2, 45), std::make_tuple(2, 59), std::make_tuple(2, 63),
        std::make_tuple(2, 71), std::make_tuple(2, 78), std::make_tuple(2, 80),
        std::make_tuple(2, 81), std::make_tuple(2, 86), std::make_tuple(2, 88),
        std::make_tuple(2, 89), std::make_tuple(2, 104), std::make_tuple(2, 109),
        std::make_tuple(2, 110), std::make_tuple(2, 114), std::make_tuple(2, 115),
        std::make_tuple(3, 14), std::make_tuple(3, 15), std::make_tuple(3, 16),
        std::make_tuple(3, 23), std::make_tuple(3, 24), std::make_tuple(3, 25),
        std::make_tuple(3, 32), std::make_tuple(3, 33), std::make_tuple(3, 34),
        std::make_tuple(3, 50), std::make_tuple(3, 54), std::make_tuple(3, 67),
        std::make_tuple(3, 70), std::make_tuple(3, 92), std::make_tuple(3, 95),
        std::make_tuple(3, 97), std::make_tuple(3, 98), std::make_tuple(3, 100),
        std::make_tuple(3, 102), std::make_tuple(3, 103), std::make_tuple(3, 118),
        std::make_tuple(3, 119), std::make_tuple(3, 121), std::make_tuple(3, 122),
        std::make_tuple(3, 124), std::make_tuple(4, 8), std::make_tuple(4, 11),
        std::make_tuple(4, 14), std::make_tuple(4, 17), std::make_tuple(4, 20),
        std::make_tuple(4, 23), std::make_tuple(4, 26), std::make_tuple(4, 29),
        std::make_tuple(4, 32), std::make_tuple(4, 37), std::make_tuple(4, 42),
        std::make_tuple(4, 48), std::make_tuple(4, 52), std::make_tuple(4, 74),
        std::make_tuple(4, 75), std::make_tuple(4, 79), std::make_tuple(4, 80),
        std::make_tuple(4, 83), std::make_tuple(4, 87), std::make_tuple(4, 88),
        std::make_tuple(4, 93), std::make_tuple(4, 94), std::make_tuple(4, 97),
        std::make_tuple(4, 99), std::make_tuple(4, 102), std::make_tuple(5, 10),
        std::make_tuple(5, 13), std::make_tuple(5, 16), std::make_tuple(5, 19),
        std::make_tuple(5, 22), std::make_tuple(5, 25), std::make_tuple(5, 28),
        std::make_tuple(5, 31), std::make_tuple(5, 34), std::make_tuple(5, 57),
        std::make_tuple(5, 61), std::make_tuple(5, 66), std::make_tuple(5, 69),
        std::make_tuple(5, 105), std::make_tuple(5, 108), std::make_tuple(5, 110),
        std::make_tuple(5, 111), std::make_tuple(5, 113), std::make_tuple(5, 115),
        std::make_tuple(5, 116), std::make_tuple(5, 117), std::make_tuple(5, 120),
        std::make_tuple(5, 121), std::make_tuple(5, 123), std::make_tuple(5, 124) };

    REQUIRE(mesh->BoundaryPoints() == all_boundaries_true);

  }

  SECTION("Mesh with at least some multi physics") {

    PetscOptionsSetValue(NULL, "--mesh-file", "fluid_layer_over_elastic_cartesian_2D_50s.e");
    PetscOptionsSetValue(NULL, "--model-file", "fluid_layer_over_elastic_cartesian_2D_50s.e");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    auto mesh = Mesh::Factory(options);
    mesh->read();
    std::unique_ptr<ExodusModel> model(new ExodusModel(options));
    model->read();
    mesh->setupTopology(model, options);
    std::unique_ptr<Problem> problem(Problem::Factory(options));
    auto elements = problem->initializeElements(mesh, model, options);
    mesh->setupGlobalDof(elements[0], options);
    std::set<std::string> all_fields_true = { "2delastic", "fluid" };
    REQUIRE(mesh->AllFields() == all_fields_true);

  }

  SECTION("Mesh independent functionality.") {

    PetscOptionsClear(NULL);
    REQUIRE(Mesh::numFieldPerPhysics("fluid") == 1);
    REQUIRE(Mesh::numFieldPerPhysics("2delastic") == 2);
    REQUIRE(Mesh::numFieldPerPhysics("3delastic") == 3);
    REQUIRE_THROWS_AS(Mesh::numFieldPerPhysics("4delastic"), std::runtime_error);

  }

}
