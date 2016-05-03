#include <iostream>
#include <Physics/AcousticNew.h>
#include <Element/Element.h>
#include <Element/HyperCube/QuadNew.h>
#include "catch.h"

TEST_CASE("test_templates", "[template]") {

  std::cout << "Hello world." << std::endl;
  auto elm = std::make_shared<QuadNew<AcousticNew>>();
  auto elm1 = elm->clone();
  elm1->prepareStiffness();

}