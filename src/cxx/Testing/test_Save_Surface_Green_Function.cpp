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

typedef SaveSurfaceGreenFunction<SaveSurfaceGreenFunctionTestStub> TestStub;

TEST_CASE("Test the saving of Green functions (for noise tomography)", "[surface_green]") {

  SECTION("Constructor") {

    PetscOptionsClear(NULL);
    PetscOptionsSetValue(NULL, "--testing", "true");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();
    REQUIRE_THROWS_AS(new TestStub(options), std::runtime_error);

    PetscOptionsSetValue(NULL, "--duration", "1");
    PetscOptionsSetValue(NULL, "--time-step", "1");
    options->setOptions();
    REQUIRE_NOTHROW(std::unique_ptr<TestStub> (new TestStub(options)));

  }

  SECTION("recordDynamicFields") {

    PetscOptionsClear(NULL);
    PetscOptionsSetValue(NULL, "--testing", "true");
    PetscOptionsSetValue(NULL, "--duration", "10");
    PetscOptionsSetValue(NULL, "--time-step", "1");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    /* Generate SaveSurface test elements. */
    std::vector<std::unique_ptr<TestStub>> elements;
    for (PetscInt i = 0; i < 10; i++) {
      elements.emplace_back(new TestStub(options));
    }

    /* Require to generate the proper sub communicator. */
    int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int size; MPI_Comm_size(PETSC_COMM_WORLD, &size);
    std::vector<PetscInt> ranks(size); ranks[rank] = rank;
    utilities::appendToGlobalRankTags(ranks, "SaveSurface");

    RealMat mock_field(TestStub::NumIntPnt(), TestStub::PullElementalFields().size());
    for (PetscInt i = 0; i < options->NumTimeSteps(); i++) {
      mock_field.col(0).setConstant(2 * i); mock_field.col(1).setConstant(2 * i + 1);
      for (auto &elm: elements) {
        elm->recordDynamicFields(mock_field);
      }
    }
  }


}
