#include <iostream>
#include <Physics/Scalar.h>
#include <Utilities/Options.h>
#include <Element/HyperCube/QuadP1.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Physics/SaveSurfaceGreenFunction.h>
#include "catch.h"

TEST_CASE("Test the saving of Green functions (for noise tomography)", "[surface_green]") {

  PetscOptionsClear(NULL);
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--polynomial-order", "4",
      NULL
  };

  /* Fake setting via command line. */
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(NULL, &argc, &argv, NULL);
  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  SaveSurfaceGreenFunction<Scalar<TensorQuad<QuadP1>>> elm(options);
}