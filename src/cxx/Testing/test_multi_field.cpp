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

  mesh->setupGlobalDof(1, 1, 1, 0, 2, model);

  PetscInt csize;
  PetscInt *idx;

  // Use this for Csize.
  DMPlexGetClosureIndices(mesh->DistributedMesh(), mesh->MeshSection(), 0, &csize, &idx);
  DMPlexGetClosureIndices(mesh->DistributedMesh(), mesh->MeshSection(), 0, &csize, &idx);
  std::vector<PetscInt> idx_0(idx, idx+csize);
  DMPlexGetClosureIndices(mesh->DistributedMesh(), mesh->MeshSection(), 1, &csize, &idx);
  std::vector<PetscInt> idx_1(idx, idx+csize);
  DMPlexGetClosureIndices(mesh->DistributedMesh(), mesh->MeshSection(), 2, &csize, &idx);
  std::vector<PetscInt> idx_2(idx, idx+csize);
  DMPlexGetClosureIndices(mesh->DistributedMesh(), mesh->MeshSection(), 3, &csize, &idx);
  std::vector<PetscInt> idx_3(idx, idx+csize);

  std::vector<double> dat_0(idx_0.size()), dat_1(idx_1.size()), dat_2(idx_2.size()), dat_3(idx_3.size());
  std::fill(dat_0.begin(), dat_0.end(), -1.0);
  std::fill(dat_1.begin(), dat_1.end(), 1.0);
  std::fill(dat_2.begin(), dat_2.end(), 2.0);
  std::fill(dat_3.begin(), dat_3.end(), 3.0);


  std::cout << idx_0.size() << ' '  << idx_2.size() << std::endl;
  Vec test_dm;
  DMCreateLocalVector(mesh->DistributedMesh(), &test_dm);
//  VecSetValues(test_dm, idx_0.size(), idx_0.data(), dat_0.data(), INSERT_VALUES);
  VecSetValues(test_dm, idx_3.size(), idx_3.data(), dat_3.data(), INSERT_VALUES);
//  VecSetValues(test_dm, idx_2.size(), idx_2.data(), dat_2.data(), INSERT_VALUES);



  std::vector<double> dat_scratch(idx_1.size());
  VecGetValues(test_dm, idx_1.size(), idx_1.data(), dat_scratch.data());
  for (auto i : dat_scratch) { std::cout << "DAT: " << i << std::endl; }
  std::cout << "SIZE: " << dat_scratch.size() << std::endl;

//  std::cout << mesh->getElementCoordinateClosure(13) << std::endl;




}