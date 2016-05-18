#include "catch.h"
#include <petsc.h>
#include <Utilities/Options.h>
#include <Mesh/Mesh.h>

extern "C" {
#include <Utilities/PETScExtensions.h>
#include <Utilities/kdtree.h>
}

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

  Mesh *mesh = Mesh::factory(options);
  mesh->read(options);

  ExodusModel *model = new ExodusModel(options);
  model->initializeParallel();

  mesh->setupGlobalDof(1, 3, 9, 0, 2, model);

  PetscScalar *buffer;
  PetscInt bsize;
  Vec test_dm;
  DMCreateLocalVector(mesh->DistributedMesh(), &test_dm);
//  DMPlexGetClosureWorkArray(mesh->DistributedMesh(), 4, mesh->MeshSection(), &bsize, &buffer);

  PetscInt csize;
  PetscInt *idx;
  DMPlexGetClosureIndices(mesh->DistributedMesh(), mesh->MeshSection(), 0, &csize, &idx);


}