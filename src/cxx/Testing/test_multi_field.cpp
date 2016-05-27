#include "catch.h"
#include <petsc.h>
#include <Utilities/Options.h>
#include <Mesh/Mesh.h>
#include <Model/ExodusModel.h>
#include <iostream>

extern "C" {
#include <Utilities/PETScExtensions.h>
#include <Utilities/kdtree.h>
}

using namespace std;

TEST_CASE("test_multi_field", "[multi_field]") {

  std::string e_file = "../../salvus_data/unit_test_meshes/fluid_layer_over_elastic_cartesian_2D_50s.e";
//  std::string e_file = "../../salvus_data/unit_test_meshes/simple_quadmesh_2x2.e";

  // Set up custom command line arguments.
  PetscOptionsClear();
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--exodus_file_name", e_file.c_str(),
      "--exodus_model_file_name", e_file.c_str(),
      "--mesh_type", "newmark",
      "--polynomial_order", "4"
  };

  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(&argc, &argv, NULL);

  Options options;
  options.setOptions();

  ExodusModel *model = new ExodusModel(options);
  model->initializeParallel();

  Mesh *mesh = Mesh::factory(options);
  mesh->read(options);
  mesh->setupGlobalDof(1, 1, 1, 0, 2, model);

  // Setup elements.
  for (PetscInt i = 0; i < mesh->NumberElementsLocal(); i++) {
    std::cout << "ELEMENT: " << i << "\n";
    for (auto f : mesh->ElementFields(i)) { std::cout << f << std::endl; }
  }

}