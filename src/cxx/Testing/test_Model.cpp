#include <iostream>
#include <salvus.h>
#include "catch.h"


TEST_CASE("Unit test model", "[model]") {

  PetscOptionsClear(NULL);
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--polynomial-order", "1", NULL};

  /* Fake setting via command line. */
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(NULL, &argc, &argv, NULL);

  SECTION("Test quads") {

    PetscOptionsSetValue(NULL, "--model-file", "quad_eigenfunction.e");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    std::unique_ptr<ExodusModel> model(new ExodusModel(options));
    model->read();

    /* Test coordinates to get first element center. */
    RealVec2 test_center; test_center << 0.0, 0.0;

    SECTION("Get correct parameters on element zero") {

      for (auto i : {0, 1, 2, 3}) {
        REQUIRE(model->getElementalMaterialParameterAtVertex(
            test_center, "VPV", 0) == Approx(5800));
        REQUIRE(model->getElementalMaterialParameterAtVertex(
            test_center, "VS", 0) == Approx(0));
      }

    }

    SECTION("Get correct element type (fluid)") {

      REQUIRE(model->getElementType(test_center) == "fluid");

    }

    SECTION("Fail because a parameter does not exist.") {

      for (auto i : {0, 1, 2, 3}) {
        REQUIRE_THROWS_AS(model->getElementalMaterialParameterAtVertex(
            test_center, "korbinian", 0), std::runtime_error);
      }
    }

    SECTION("Fail because no nodal variables are defined.") {
      /* Only elemental variables defined for this mesh. */
      REQUIRE_THROWS_AS(model->getNodalParameterAtNode(test_center, "fail"),
                        std::runtime_error);
    }

    SECTION("Check side sets") {

      std::vector<std::string> true_side_sets { "x0", "x1", "y0", "y1" };
      for (auto i: {0, 1, 2, 3}) {
        REQUIRE(model->SideSetName(i) == true_side_sets[i]);
      }
      REQUIRE_THROWS_AS(model->SideSetName(4), std::runtime_error);

    }

  }

  SECTION("Test hexes") {

    PetscOptionsSetValue(NULL, "--model-file", "hex_eigenfunction.e");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    std::unique_ptr<ExodusModel> model(new ExodusModel(options));
    model->read();

    /* Test coordinates to get first element center. */
    RealVec3 test_center; test_center << 0.0, 0.0, 0.0;

    SECTION("Get correct parameters on element zero") {

      for (auto i : {0, 1, 2, 3, 4, 5, 6, 7}) {
        REQUIRE(model->getElementalMaterialParameterAtVertex(
            test_center, "VPV", 0) == Approx(5800));
        REQUIRE(model->getElementalMaterialParameterAtVertex(
            test_center, "VS", 0) == Approx(0));
      }

    }

    SECTION("Get correct element type (fluid)") {

      REQUIRE(model->getElementType(test_center) == "fluid");

    }

    SECTION("Fail because a parameter does not exist.") {

      for (auto i : {0, 1, 2, 3, 4, 5, 6, 7}) {
        REQUIRE_THROWS_AS(model->getElementalMaterialParameterAtVertex(
            test_center, "korbinian", 0), std::runtime_error);
      }
    }

    SECTION("Fail because no nodal variables are defined.") {
      /* Only elemental variables defined for this mesh. */
      REQUIRE_THROWS_AS(model->getNodalParameterAtNode(test_center, "fail"),
                        std::runtime_error);
    }

    SECTION("Check side sets") {

      std::vector<std::string> true_side_sets{"x0", "x1", "y0", "y1", "z0", "z1"};
      for (auto i: {0, 1, 2, 3, 4, 5}) {
        REQUIRE(model->SideSetName(i) == true_side_sets[i]);
      }
      REQUIRE_THROWS_AS(model->SideSetName(6), std::runtime_error);

    }

  }

  SECTION("Test nodal pars") {

    PetscOptionsSetValue(NULL, "--model-file", "nodal_hex.e");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    std::unique_ptr<ExodusModel> model(new ExodusModel(options));
    model->read();

    SECTION("Correct values") {

      /* Test coordinates to get first element center. */
      RealVec3 test_center;

      /* Node VP values were written as node numbers themselves. */
      test_center << 0.0, 0.0, 0.0;
      REQUIRE(model->getNodalParameterAtNode(test_center, "VP") == Approx(0));

      test_center << 0.0, 0.0, 0.5;
      REQUIRE(model->getNodalParameterAtNode(test_center, "VP") == Approx(1));

      test_center << 0.0, 0.5, 0.5;
      REQUIRE(model->getNodalParameterAtNode(test_center, "VP") == Approx(4));

      test_center << 0.0, 0.5, 0.0;
      REQUIRE(model->getNodalParameterAtNode(test_center, "VP") == Approx(3));

      test_center << 0.5, 0.0, 0.0;
      REQUIRE(model->getNodalParameterAtNode(test_center, "VP") == Approx(9));

      test_center << 0.5, 0.0, 0.5;
      REQUIRE(model->getNodalParameterAtNode(test_center, "VP") == Approx(10));

      test_center << 0.5, 0.5, 0.5;
      REQUIRE(model->getNodalParameterAtNode(test_center, "VP") == Approx(13));

      test_center << 0.5, 0.5, 0.0;
      REQUIRE(model->getNodalParameterAtNode(test_center, "VP") == Approx(12));

      /* Only VP should have correct values. */
      test_center << 0.0, 0.0, 0.5;
      REQUIRE(model->getNodalParameterAtNode(test_center, "RHO") == Approx(0));

    }

    SECTION("Throw if trying to get nodal param that doens't exist") {

      RealVec3 test_center;
      test_center << 0.0, 0.0, 0.0;

      REQUIRE_THROWS_AS(model->getNodalParameterAtNode(
          test_center, "VPP"), std::runtime_error);

    }

    SECTION("Throw if trying to get elemental params") {

      RealVec3 test_center;
      test_center << 0.0, 0.0, 0.0;

      REQUIRE_THROWS_AS(model->getElementalMaterialParameterAtVertex(
          test_center, "VP", 0), std::runtime_error);

    }

  }

  SECTION("Test read/write") {

    PetscOptionsSetValue(NULL, "--model-file", "quad_eigenfunction.e");
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    std::unique_ptr<ExodusModel> model(new ExodusModel());
    SECTION("Fail if Exodus file name is not specified") {    
      REQUIRE_THROWS_AS(model->read(), std::runtime_error);
    }

    model->setExodusFilename( options->ModelFile() );
    model->read();
    

    SECTION("Write exodus model to disk") {
      std::string output_filename("moep.e");
      model->write(output_filename.c_str());
      std::ifstream f(output_filename.c_str());
      REQUIRE(f.good());
    }
  }

}
