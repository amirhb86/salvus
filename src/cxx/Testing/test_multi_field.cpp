#include "catch.h"
#include <petsc.h>
#include <Utilities/Options.h>
#include <Mesh/Mesh.h>
#include <Model/ExodusModel.h>

TEST_CASE("test_multi_field", "[multi_field]") {

  std::string e_file = "../../salvus_data/unit_test_meshes/fluid_layer_over_elastic_cartesian_2D_50s.e";

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

  Mesh *mesh = Mesh::factory(options);
  mesh->read(options);

  ExodusModel *model = new ExodusModel(options);
  model->initializeParallel();

  Eigen::VectorXd pnt;
  pnt.resize(2); pnt(0) = 1e5; pnt(1) = 1e5;
  std::cout << model->getElementType(pnt) << std::endl;

}