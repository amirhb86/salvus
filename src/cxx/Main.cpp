#include <iostream>
#include <vector>
#include <mpi.h>
#include <petscsys.h>
#include <memory>

#include "Mesh/Mesh.h"
#include "Problem/Problem.h"

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
  std::shared_ptr<Element> reference_element = Element::factory(options);

  // Use above elements to define the problem.
  Problem *problem = Problem::factory(options.ProblemType());
  problem->initialize(mesh, model, reference_element, options);
  problem->solve(options);

  PetscFinalize();
}
