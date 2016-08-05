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
    std::unique_ptr<Options> options(new Options);
    options->setOptions();
    auto mesh = Mesh::Factory(options);
    REQUIRE_THROWS_AS(mesh->read(), std::runtime_error);
  }

  SECTION("Correctly initialize DM (2D quad).") {

    PetscOptionsSetValue(NULL, "--mesh-file", "quad_eigenfunction.e");
    PetscOptionsSetValue(NULL, "--model-file", "quad_eigenfunction.e");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    /* Ensure mesh was read. */
    auto mesh = Mesh::Factory(options);
    mesh->read();
    REQUIRE(mesh->DistributedMesh() != NULL);

    /* Attach model and global dofs. */
    std::unique_ptr<ExodusModel> model(new ExodusModel(options));
    model->read();
    mesh->setupGlobalDof(model, options);
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
    REQUIRE(mesh->TotalCouplingFields(0)[0] == "boundary");

    /* Ensure that if we're in the middle of a homogeneous section, there's no coupling. */
    REQUIRE(mesh->TotalCouplingFields(5).empty());

    /* Make sure we know we're a quad. */
    REQUIRE(mesh->baseElementType() == "quad");

    /* Make sure we have the correct number of elements. */
    REQUIRE(mesh->NumberElementsLocal() == 16);

    /* Ensure that we get a correct listing of all the boundary points. */
    std::set<PetscInt> all_boundaries_true {
        41, 44, 45, 48, 51, 52, 56, 61, 65, 70, 73, 74, 76, 78, 79, 80 };
    REQUIRE(mesh->BoundaryPoints() == all_boundaries_true);



  }

  SECTION("Correctly initialize DM (3D hex).") {
    PetscOptionsSetValue(NULL, "--mesh-file",  "small_hex_mesh_to_test_sources.e");
    PetscOptionsSetValue(NULL, "--model-file", "small_hex_mesh_to_test_sources.e");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    /* Ensure mesh was read. */
    auto mesh = Mesh::Factory(options);
    mesh->read();
    REQUIRE(mesh->DistributedMesh() != NULL);

    /* Attach model and global dofs. */
    std::unique_ptr<ExodusModel> model(new ExodusModel(options));
    model->read();
    mesh->setupGlobalDof(model, options);
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
    REQUIRE(mesh->TotalCouplingFields(0)[0] == "boundary");
    REQUIRE(mesh->TotalCouplingFields(1)[0] == "boundary");

    /* TODO: See above TODO. */
    /* Ensure that if we're in the middle of a homogeneous section, there's no coupling. */
    /* REQUIRE(mesh->TotalCouplingFields(5).empty() == true); */

    /* Ensure we know we're a hex. */
    REQUIRE(mesh->baseElementType() == "hex");

    /* Ensure we have the proper number of elements. */
    REQUIRE(mesh->NumberElementsLocal() == 8);

    /* Ensure that we get a correct listing of all the boundary points. */
    std::set<PetscInt> all_boundaries_true {
        35, 37, 40, 41, 42, 45, 46, 48, 50, 51, 52, 54, 55, 57, 59, 60,
        61, 63, 64, 66, 67, 68, 69, 70 };
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
    mesh->setupGlobalDof(model, options);
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
