#include <iostream>
#include <Physics/Scalar.h>
#include <Utilities/Utilities.h>
#include <Utilities/Options.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Physics/SaveSurfaceGreenFunction.h>
#include "catch.h"
#include <hdf5_hl.h>

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

    /* Require to generate the proper sub communicator. */
    int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int size; MPI_Comm_size(PETSC_COMM_WORLD, &size);

    /* Generate SaveSurface test elements. */
    std::vector<std::unique_ptr<TestStub>> elements;
    for (PetscInt i = 0; i < 1 * (rank + 1); i++) {
      if (rank != 1) { elements.emplace_back(new TestStub(options)); }
    }

    /* Make a mixin on all ranks but 1. */
    PetscInt rsize = size > 1 ? size - 1 : size;
    std::vector<PetscInt> ranks(rsize);
    for (PetscInt i = 0, j = 0; i < size; i++) { if (i != 1) { ranks[j++] = i; }}
    utilities::appendToGlobalRankTags(ranks, "SaveSurface");

    /* Pretend we have some funny time loop. */
    RealMat mock_field(TestStub::NumIntPnt(), TestStub::PullElementalFields().size());
    for (PetscInt i = 0; i < options->NumTimeSteps(); i++) {
      mock_field.col(0).setConstant(rank);
      mock_field.col(1).setConstant(2 * i + 1);
      for (auto &elm: elements) {
        elm->recordDynamicFields(mock_field);
      }
    }

  }

  SECTION("check if write worked.") {

    PetscOptionsClear(NULL);
    PetscOptionsSetValue(NULL, "--testing", "true");
    PetscOptionsSetValue(NULL, "--duration", "10");
    PetscOptionsSetValue(NULL, "--time-step", "1");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    /* Require to generate the proper sub communicator. */
    int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int size; MPI_Comm_size(PETSC_COMM_WORLD, &size);

    PetscInt tot_elem = 0;
    for (PetscInt i = 0; i < 1 * (rank + 1); i++) { if (rank != 1) { tot_elem++; } }
    MPI_Allreduce(MPI_IN_PLACE, &tot_elem, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD);

    /* Open to see if write worked. */
    PetscInt num_c = SaveSurfaceGreenFunctionTestStub::PullElementalFields().size();
    PetscInt num_p = SaveSurfaceGreenFunctionTestStub::NumIntPnt();
    UncompressedWavefieldContainer<PetscReal> t(
        options->NumTimeSteps(), tot_elem,
        SaveSurfaceGreenFunctionTestStub::PullElementalFields().size(),
        SaveSurfaceGreenFunctionTestStub::NumIntPnt());

    /* Can we get some sort of MD5 hash here instead? */
    hid_t fileid = H5Fopen("test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    H5LTread_dataset(fileid, "wavefield_data", H5T_NATIVE_DOUBLE, &t.data());
    if (size == 4) {
      REQUIRE(t(8, 0, 0, 0) == Approx(0));
      REQUIRE(t(8, 1, 0, 0) == Approx(2));
      REQUIRE(t(8, 4, 0, 0) == Approx(3));
    } else if (size == 2) {
      REQUIRE(t(4, 0, 1, 0) == Approx(9));
      REQUIRE(t(6, 1, 0, 0) == Approx(0));
      REQUIRE(t(3, 4, 1, 0) == Approx(15));
    } else if (size == 1) {
      REQUIRE(t(4, 0, 1, 0) == Approx(9));
      REQUIRE(t(6, 1, 1, 0) == Approx(15));
      REQUIRE(t(3, 4, 1, 0) == Approx(15));
    }


  }


}
