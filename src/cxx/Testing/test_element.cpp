#include <iostream>
#include <memory>

#include "catch.h"
#include <petsc.h>

#include <Element/Element.h>
#include <Element/ElementAdapter.h>
#include <Element/HyperCube/Quad.h>
#include <Element/HyperCube/QuadP1.h>
#include <Physics/Elastic2D.h>
#include <Physics/Acoustic2D.h>
#include <Utilities/Options.h>
#include <Element/ElementAdapter.h>

TEST_CASE("element", "[element]") {

  typedef class ElementAdapter<Acoustic2D<Quad<QuadP1>>> AcousticQuadP1;

  PetscOptionsClear();
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--element_shape", "quad_new", NULL};

  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(&argc, &argv, NULL);

  std::unique_ptr<Options> options;
  options->setOptions();

  auto elm = Element::Factory({"u"}, {}, options);

  /**** CHECK CLASS NAMES HERE!!! ****/

}