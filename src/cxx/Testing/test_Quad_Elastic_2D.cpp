#include <iostream>
#include <Utilities/Types.h>
#include <Mesh/Mesh.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Physics/Scalar.h>
#include <Element/HyperCube/QuadP1.h>
#include <Element/ElementAdapter.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Problem/Problem.h>
#include <petscviewerhdf5.h>
#include "catch.h"

TEST_CASE("Test elastic in 2D", "[elastic_2d]") {

  std::string e_file = "elastic_test.e";

  PetscOptionsClear(NULL);
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--mesh-file", e_file.c_str(),
      "--model-file", e_file.c_str(),
      "--polynomial-order", "4",
      "--time-step", "1e-2",
      "--duration", "0.1",
      "--number-of-sources", "2",
      "--source-type", "ricker",
      "--source-location-x", "50000,90000",
      "--source-location-y", "50000,90000",
      "--ricker-amplitude", "100,100",
      "--ricker-time-delay", "1.0,1.5",
      "--ricker-center-freq", "0.5,0.5",
      NULL };

  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(NULL, &argc, &argv, NULL);

  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  std::unique_ptr<Problem>      problem(Problem::Factory(options));
  std::unique_ptr<ExodusModel>  model(new ExodusModel(options));
  std::unique_ptr<Mesh>         mesh(Mesh::Factory(options));

  model->read();
  mesh->read();

  /* Setup topology from model and mesh. */
  mesh->setupTopology(model, options);

  /* Setup elements from model and topology. */
  auto elements = problem->initializeElements(mesh, model, options);

  /* Setup global degrees of freedom based on element 0. */
  mesh->setupGlobalDof(elements[0], options);

  auto fields = problem->initializeGlobalDofs(elements, mesh);

  PetscReal time = 0;
  PetscInt time_idx = 0;
  while (time < options->Duration()) {

    std::tie(elements, fields) = problem->assembleIntoGlobalDof(
        std::move(elements), std::move(fields), time, time_idx,
        mesh->DistributedMesh(), mesh->MeshSection(), options);

    fields = problem->applyInverseMassMatrix(std::move(fields));
    std::tie(fields, time) = problem->takeTimeStep
        (std::move(fields), time, options);

    time_idx++;
    std::cout << "TIME:      " << time << "\r"; std::cout.flush();

  }

  /* If the max value, and global index, is the same, assume that the regression has passed. */
//  PetscInt ind; PetscReal max;
//  PetscInt regression_ind = 17744; PetscReal regression_max = 2.77143e-07;
//  VecMax(fields["u"]->mGlb, &ind, &max);
//  REQUIRE(max == Approx(regression_max));
//  REQUIRE(ind == regression_ind);

}
