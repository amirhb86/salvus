#include "catch.h"
#include <petsc.h>
#include <Utilities/Options.h>
#include <Problem/ProblemNew.h>
#include <Model/Model.h>
#include <Mesh/ElasticAcousticNewmark3D.h>

TEST_CASE("Test new problem formulation", "[problem_new]") {

  PetscOptionsClear();
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--element_shape", "quad_new",
      "--polynomial_order", "4", NULL
  };

  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(&argc, &argv, NULL);

  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  std::unique_ptr<Mesh> mesh(new ElasticAcousticNewmark3D());
  std::unique_ptr<Model> model(new Model());
  std::unique_ptr<ProblemNew> problem(ProblemNew::Factory(options));

  auto elements = problem->initializeElements(mesh, model, options);


}