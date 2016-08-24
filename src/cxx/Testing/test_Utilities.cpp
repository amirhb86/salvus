#include "catch.h"
#include <salvus.h>
#include <Utilities/Utilities.h>
#include <Physics/SaveSurfaceGreenFunction.h>

TEST_CASE("Unit test utilities", "[utilities]") {

  SECTION("getRanksForMixin") {

    PetscOptionsClear(NULL);
    PetscOptionsSetValue(NULL, "--testing", "true");
    PetscOptionsSetValue(NULL, "--polynomial-order", "4");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    std::vector<std::unique_ptr<Element>> elements;
    int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int size; MPI_Comm_size(PETSC_COMM_WORLD, &size);

    if (rank != 1) {
      elements.emplace_back(
          new ElementAdapter<
              SaveSurfaceGreenFunction<
                  Scalar<
                      TensorQuad<
                          QuadP1>>>>( options));
    } else {
      elements.emplace_back(
          new ElementAdapter<
              Scalar<
                  TensorQuad<
                      QuadP1>>>(options));
    }

    /* If we're running in serial. */
    if (size == 1) {
      std::vector<PetscInt> serial {0};
      auto ranks = utilities::getRanksForMixin(elements, "SaveSurface");
      utilities::appendToGlobalRankTags(ranks, "SaveSurface");
      REQUIRE(ranks == serial);
      REQUIRE(utilities::GetWorldRanksForTag("SaveSurface") == serial);
    } else {
      std::vector<PetscInt> parallel(size - 1);
      for (PetscInt i = 0, j = 0; j < parallel.size(); i++) { if (i != 1) { parallel[j++] = i; } }
      auto ranks = utilities::getRanksForMixin(elements, "SaveSurface");
      utilities::appendToGlobalRankTags(ranks, "SaveSurface");
      REQUIRE(ranks == parallel);
      REQUIRE(utilities::GetWorldRanksForTag("SaveSurface") == parallel);
    }

    REQUIRE_THROWS_AS(utilities::GetWorldRanksForTag("error"), std::runtime_error);

  }

  SECTION("Test string has exception") {
    REQUIRE(utilities::stringHasExtension("my_test.h5", "h5"));
    REQUIRE_FALSE(utilities::stringHasExtension("my_test.h5", "e"));
  }

  SECTION("mpitype") {
    REQUIRE(utilities::mpitype<double>() == MPI_DOUBLE);
    REQUIRE(utilities::mpitype<float>() == MPI_FLOAT);
    REQUIRE(utilities::mpitype<PetscInt>() == MPIU_INT);
  }

  SECTION("broadcastNumberFromRank") {
    REQUIRE(utilities::broadcastNumberFromRank(1.0, 0) == 1.0);
  }

  SECTION("broadcastNumberVecFromRank") {
    std::vector<float> bcast{1.0, 2.0};
    REQUIRE(utilities::broadcastNumberVecFromRank(bcast, 0) == bcast);
  }

  SECTION("broadcastStringFromRank") {
    std::string bcast("test");
    REQUIRE(utilities::broadcastStringFromRank(bcast, 0) == bcast);
  }

  SECTION("broadcastStringVecFromRank") {
    std::vector<std::string> bcast {"test1", "test2"};
    REQUIRE(utilities::broadcastStringVecFromRank(bcast, 0) == bcast);
  }

}