#define CATCH_CONFIG_RUNNER
#include "catch.h"
#include <Eigen/Dense>
#include <petsc.h>
#include "../../../include/Element/ElementAdapter.h"
#include <Element/HyperCube/Quad.h>
#include <Element/HyperCube/QuadP1.h>
#include <Physics/Acoustic2D.h>
#include <Element/Simplex/Triangle.h>
#include <Element/Simplex/TriP1.h>
#include <Element/HyperCube/Quad.h>
#include <Model/ExodusModel.h>

#include <Utilities/Logging.h>

INIT_LOGGING_STATE();

int main(int argc, char *argv[]) {

  // Init Salvus command line arguments.
  PetscInitialize(&argc, &argv, NULL, NULL);

  // Run all unit tests.
  int result = Catch::Session().run(argc, argv);

  // Clean up PETSc.
  PetscFinalize();

  return result;
}

