#include "catch.h"
#include <Element/Element.h>
#include <Element/HyperCube/Quad/Acoustic.h>
#include <Eigen/Dense>

Quad *setup_simple_quad_source(Options options) {

  // Simple model.
  ExodusModel *model = new ExodusModel(options);
  model->initializeParallel();

  // Get element from options.
  Element *reference_element = Element::factory(options);
  Quad *reference_quad = dynamic_cast<Quad*> (reference_element);

  // Make things easy by assuming a reference element.
  // NOTE THE ELEMENT IS DISTORTED x -> [-2, 1], y -> [-6, 1]
  Eigen::Matrix<double,2,4> coord;
  coord << -2, +1, +1, -2,
      -6, -6, +1, +1;
  reference_quad->SetVtxCrd(coord);
  reference_quad->attachMaterialProperties(model);
  reference_quad->setupGradientOperator();

  return reference_quad;

}

TEST_CASE("Test whether source interpolates properly.", "[source]") {

  Options options;
  options.__SetPolynomialOrder(3);
  options.__SetExodusMeshFile("simple_quadmesh_2x2.e");
  options.__SetExodusModelFile("simple_quadmesh_2x2.e");
  options.__SetElementShape("quad");
  options.__SetPhysicsSystem("acoustic");
//  options.__

  auto sources = Source::factory(options);

}