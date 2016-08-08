#include <iostream>
#include <salvus.h>
#include "catch.h"


TEST_CASE("Test to make sure proper elements are returned.", "[element]") {

  PetscOptionsClear(NULL);
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--polynomial-order", "1", NULL };

  /* Fake setting via command line. */
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(NULL, &argc, &argv, NULL);
  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  REQUIRE(Element::Factory("quad", {"fluid"}, {}, options)->Name() ==
      "Scalar_TensorQuad_QuadP1");
  REQUIRE(Element::Factory("quad", {"fluid"}, {"boundary"}, options)->Name() ==
      "HomogeneousDirichlet_Scalar_TensorQuad_QuadP1");
  REQUIRE(Element::Factory("quad", {"fluid"}, {"2delastic"}, options)->Name() ==
      "SolidToFluid2D_Scalar_TensorQuad_QuadP1");
  REQUIRE(Element::Factory("quad", {"fluid"}, {"2delastic", "boundary"}, options)->Name() ==
      "HomogeneousDirichlet_SolidToFluid2D_Scalar_TensorQuad_QuadP1");
  REQUIRE(Element::Factory("quad", {"2delastic"}, {}, options)->Name() ==
      "Elastic2D_TensorQuad_QuadP1");
  REQUIRE(Element::Factory("quad", {"2delastic"}, {"fluid"}, options)->Name() ==
      "FluidToSolid2D_Elastic2D_TensorQuad_QuadP1");
  REQUIRE(Element::Factory("quad", {"2delastic"}, {"boundary"}, options)->Name() ==
      "HomogeneousDirichlet_Elastic2D_TensorQuad_QuadP1");
  REQUIRE(Element::Factory("quad", {"2delastic"}, {"fluid", "boundary"}, options)->Name() ==
      "HomogeneousDirichlet_FluidToSolid2D_Elastic2D_TensorQuad_QuadP1");

  REQUIRE(Element::Factory("hex", {"fluid"}, {}, options)->Name() ==
      "Scalar_TensorHex_HexP1");
  REQUIRE(Element::Factory("hex", {"fluid"}, {"boundary"}, options)->Name() ==
      "HomogeneousDirichlet_Scalar_TensorHex_HexP1");
  REQUIRE(Element::Factory("hex", {"3delastic"}, {}, options)->Name() ==
      "Elastic3D_TensorHex_HexP1");
  REQUIRE(Element::Factory("hex", {"3delastic"}, {"boundary"}, options)->Name() ==
      "HomogeneousDirichlet_Elastic3D_TensorHex_HexP1");

  /* Make sure dumb values are not allowed. */
  REQUIRE_THROWS_AS(Element::Factory("hex", {"2delastic"}, {}, options)->Name(),
                    std::runtime_error);
  REQUIRE_THROWS_AS(Element::Factory("quad", {"3delastic"}, {}, options)->Name(),
                    std::runtime_error);

  REQUIRE_THROWS_AS(Element::Factory("qua", {"2delastic"}, {}, options)->Name(),
                    std::runtime_error);
  REQUIRE_THROWS_AS(Element::Factory("quad", {"2deastic"}, {}, options)->Name(),
                    std::runtime_error);

  REQUIRE_THROWS_AS(Element::Factory("quad", {"2delastic", "fluid", "boundary"}, {}, options)->Name(),
                    std::runtime_error);
  REQUIRE_THROWS_AS(Element::Factory("hex", {"3delastic", "fluid", "boundary"}, {}, options)->Name(),
                    std::runtime_error);

  REQUIRE_THROWS_AS(Element::Factory("quad", {"2delastic"}, {"cat", "dog", "korbinian"}, options)->Name(),
                    std::runtime_error);
  REQUIRE_THROWS_AS(Element::Factory("hex", {"3delastic"}, {"fluid", "boundary", "dog"}, options)->Name(),
                    std::runtime_error);

}

