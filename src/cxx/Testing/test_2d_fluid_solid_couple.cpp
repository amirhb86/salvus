#include "catch.h"

#include <iostream>
#include <Eigen/Dense>

#include <petsc.h>
#include <Mesh/Mesh.h>
#include <Element/Element.h>
#include <Element/ElementAdapter.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>

using namespace Eigen;
using namespace std;

TEST_CASE("test_fluid_solid_couple", "[couple/fluid_solid]") {

  PetscOptionsClear(NULL);
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--exodus_file_name", "fluid_layer_over_elastic_cartesian_2D_50s.e",
      "--exodus_model_file_name", "fluid_layer_over_elastic_cartesian_2D_50s.e",
      "--mesh_type", "newmark",
      "--element_shape", "quad_new",
      "--polynomial_order", "4", NULL};

  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;

  PetscOptionsInsert(NULL, &argc, &argv, NULL);

  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  std::unique_ptr<Mesh> const &mesh = Mesh::Factory(options);
  mesh->read(options);

  std::unique_ptr<ExodusModel> model(new ExodusModel(options));
  model->initializeParallel();
  mesh->setupGlobalDof(2, model, options);

  std::vector<std::shared_ptr<Element>> elms;
  for (PetscInt i = 0; i < mesh->NumberElementsLocal(); i++) {
    elms.push_back(Element::Factory(mesh->ElementFields(i),
                                    mesh->TotalCouplingFields(i),
                                    options));
    elms.back()->SetNum(i);
  }

  for (auto &e: elms) {
    e->attachVertexCoordinates(mesh);
    e->setBoundaryConditions(mesh);
    e->attachMaterialProperties(model);
  }

  // Test against some analytic solution assuming constant velocity along the interface.
//  Eigen::MatrixXd vel = Eigen::MatrixXd::Constant(25, 3, 1.0);
//  REQUIRE(elms[0]->computeSurfaceIntegral(vel).sum() == Approx(-2 * 1.3e8));
//
//  vel.col(2).setConstant(2);
//  REQUIRE(elms[2]->computeSurfaceIntegral(vel).sum() == Approx(2 * 50000));

}
