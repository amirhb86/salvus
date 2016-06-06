// stl.
#include <iostream>
#include <vector>
#include <memory>

// 3rd party.
#include <mpi.h>
#include <petsc.h>
#include <Element/Element.h>

#include <Mesh/Mesh.h>
#include <Source/Source.h>
#include <Problem/Problem.h>
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>
#include <Utilities/Logging.h>

INIT_LOGGING_STATE();

static constexpr char help[] = "Welcome to Salvus.";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc, &argv, NULL, help);

  // Get command line options.
  Options options;
  options.setOptions();

  // Get mesh.
  Mesh *mesh = Mesh::factory(options);
  mesh->read(options);

  // Get model.
  ExodusModel *model = new ExodusModel(options);
  model->initializeParallel();

  // Get sources.
  std::vector<std::shared_ptr<Source>> sources = Source::factory(options);

  // Setup reference element.
  std::shared_ptr<Element> reference_element = Element::Factory(options);

  // Use above elements to define the problem.
  Problem *problem = Problem::factory(options.ProblemType());
  problem->initialize(mesh, model, reference_element, options);
  problem->solve(options);

  PetscFinalize();
}
