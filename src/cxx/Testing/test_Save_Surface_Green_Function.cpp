#include <iostream>
#include <Physics/Scalar.h>
#include <Utilities/Utilities.h>
#include <Utilities/Options.h>
#include <Element/Element.h>
#include <Element/ElementAdapter.h>
#include <Element/HyperCube/QuadP1.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Physics/SaveSurfaceGreenFunction.h>
#include "catch.h"

TEST_CASE("Test the saving of Green functions (for noise tomography)", "[surface_green]") {

  PetscOptionsClear(NULL);
  const char *arg[] = {"salvus_test", "--testing", "true", "--polynomial-order", "4", NULL};

  /* Fake setting via command line. */
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(NULL, &argc, &argv, NULL);
  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  std::vector<std::unique_ptr<Element>> elements;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank != 1) {
    elements.emplace_back(new ElementAdapter<SaveSurfaceGreenFunction<Scalar<TensorQuad<QuadP1>>>>(
        options));
  } else {
    elements.emplace_back(new ElementAdapter<Scalar<TensorQuad<QuadP1>>>(options));
  }

  auto ranks = utilities::getRanksForMixin(elements, "SaveSurface");
  utilities::appendToGlobalRankTags(ranks, "SaveSurface");

  RealMat f(1, 0);
  for (auto &elm: elements) {
    elm->recordDynamicFields(f);
  }

}

