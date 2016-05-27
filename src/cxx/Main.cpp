#include <iostream>
#include <vector>
#include <mpi.h>
#include <petscsys.h>
#include <memory>
#include "../../include/Element/Element.h"

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

//  std::cout << "HERE" << std::endl;
//  MPI_Abort(PETSC_COMM_WORLD, -1);
//  exit(0);
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
