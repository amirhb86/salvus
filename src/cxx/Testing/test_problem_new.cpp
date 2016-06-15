#include "catch.h"
#include <petsc.h>
#include <Utilities/Options.h>
#include <Problem/ProblemNew.h>
#include <Model/ExodusModel.h>
#include <Mesh/ElasticAcousticNewmark3D.h>

TEST_CASE("Test new problem formulation", "[problem_new]") {

  std::string e_file = "fluid_layer_over_elastic_cartesian_2D_50s.e";

  PetscOptionsClear();
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--element_shape", "quad_new",
      "--polynomial_order", "4",
      "--exodus_file_name", e_file.c_str(),
      "--exodus_model_file_name", e_file.c_str(), NULL
  };

  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(&argc, &argv, NULL);

  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  std::unique_ptr<ExodusModel> model(new ExodusModel(options));
  model->initializeParallel();

  std::unique_ptr<Mesh> mesh(new ElasticAcousticNewmark3D(options));
  mesh->read(options);
  mesh->setupGlobalDof(1, 3, 9, 0, 2, model);

  std::unique_ptr<ProblemNew> problem(ProblemNew::Factory(options));

  auto elements = problem->initializeElements(mesh, model, options);
  for (auto &e: elements) {
    std::cout << e->Num() << std::endl;
  }


}