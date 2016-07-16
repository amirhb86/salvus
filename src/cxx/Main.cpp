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

//  // Get command line options.
//  std::unique_ptr<Options> options;
//  options->setOptions();
//
//  // Get mesh.
//  auto mesh = Mesh::Factory(options);
//  mesh->read(options);
//
//  // Get model.
//  std::unique_ptr<ExodusModel> model(new ExodusModel(options));
//  model->initializeParallel();
//
//  // Get sources.
//  auto sources = Source::Factory(options);
//
//  // Use above elements to define the problem.
//  Problem *problem = Problem::Factory(options->ProblemType());
//
//  problem->initialize(std::move(mesh), model, options);
////  problem->solve(options);

  PetscFinalize();
}
