#include <iostream>
#include <Utilities/Types.h>
#include <Mesh/Mesh.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Physics/Scalar.h>
#include <Element/ElementAdapter.h>
#include <Problem/Problem.h>
#include <petscviewerhdf5.h>
#include <Element/HyperCube/Hexahedra.h>
#include <Element/HyperCube/HexP1.h>
#include "catch.h"

using namespace std;

template <typename Element>
class TestPlugin: public Element {

 public:
  TestPlugin<Element>(std::unique_ptr<Options> const &options): Element(options) {};

  void setupEigenfunctionTest(std::unique_ptr<Mesh> const &mesh,
                              std::unique_ptr<Options> const &options,
                              std::unique_ptr<Problem> &problem,
                              FieldDict &fields) {

    /* This is hardcoded for the unit test mesh. */
    PetscScalar x0 = 5e4, y0 = 5e4, z0 = 5e4, L = 1e5;
    RealVec pts_x, pts_y, pts_z;
    std::tie(pts_x, pts_y, pts_z) = Element::buildNodalPoints();
    RealVec un =
        (M_PI / L * (pts_x.array() - (x0 + L / 2))).sin() *
        (M_PI / L * (pts_y.array() - (y0 + L / 2))).sin() *
        (M_PI / L * (pts_z.array() - (z0 + L / 2))).sin();
    RealVec vn = RealVec::Zero(pts_x.size());
    RealVec an = RealVec::Zero(pts_x.size());
    problem->insertElementalFieldIntoMesh("u", Element::ElmNum(), Element::ClsMap(), un,
                                          mesh->DistributedMesh(), mesh->MeshSection(),
                                          fields);
    problem->insertElementalFieldIntoMesh("v", Element::ElmNum(), Element::ClsMap(), vn,
                                          mesh->DistributedMesh(), mesh->MeshSection(),
                                          fields);
    problem->insertElementalFieldIntoMesh("a", Element::ElmNum(), Element::ClsMap(), an,
                                          mesh->DistributedMesh(), mesh->MeshSection(),
                                          fields);

  }

  PetscReal checkEigenfunctionTestNew(std::unique_ptr<Mesh> const &mesh,
                                      std::unique_ptr<Options> const &options,
                                      const PetscScalar time,
                                      std::unique_ptr<Problem> &problem,
                                      FieldDict &fields) {

    PetscScalar x0 = 5e4, y0 = 5e4, z0 = 5e4, L = 1e5;
    RealVec pts_x, pts_y, pts_z;
    std::tie(pts_x, pts_y, pts_z) = Element::buildNodalPoints();
    RealVec un_xyz =
        (M_PI / L * (pts_x.array() - (x0 + L / 2))).sin() *
        (M_PI / L * (pts_y.array() - (y0 + L / 2))).sin() *
        (M_PI / L * (pts_z.array() - (z0 + L / 2))).sin();
    PetscScalar vp = Element::ParAtIntPts("VP").mean();
    PetscScalar un_t = cos(M_PI / L * sqrt(3) * time * vp);
    RealVec exact = un_t * un_xyz;

    RealVec u = problem->getFieldOnElement(
        "u", Element::ElmNum(), Element::ClsMap(),
        mesh->DistributedMesh(), mesh->MeshSection(), fields);

    PetscScalar element_error = (exact - u).array().abs().maxCoeff();
    return element_error;

  }
};

typedef ElementAdapter<TestPlugin<Scalar<Hexahedra<HexP1>>>> test_insert;
typedef TestPlugin<Scalar<Hexahedra<HexP1>>> test_init;
typedef ElementAdapter<Scalar<Hexahedra<HexP1>>> unguard;
typedef Scalar<Hexahedra<HexP1>> raw;

TEST_CASE("Test analytic eigenfunction solution for scalar "
              "equation in 3d", "[hex_eigenfunction]") {

  std::string e_file = "hex_eigenfunction.e";

  PetscOptionsClear(NULL);
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--mesh-file", e_file.c_str(),
      "--model-file", e_file.c_str(),
      "--time-step", "1e-2",
      "--polynomial-order", "2", NULL};
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(NULL, &argc, &argv, NULL);

  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  std::unique_ptr<Problem> problem(Problem::Factory(options));
  std::unique_ptr<ExodusModel> model(new ExodusModel(options));
  std::unique_ptr<Mesh> mesh(Mesh::Factory(options));

  model->read();
  mesh->read();
  mesh->setupGlobalDof(model, options);

  std::vector<std::unique_ptr<Element>> test_elements;
  auto elements = problem->initializeElements(mesh, model, options);
  auto fields = problem->initializeGlobalDofs(elements, mesh);

  /* Rip apart elements and insert testing mixin. */
  for (auto &e: elements) {

    /* Rip out the master Element class. */
    auto l1 = static_cast<unguard*>(e.release());

    /* Rip out the Element adapter. */
    auto l2 = static_cast<raw*>(l1);

    /* Attach the tester. */
    auto l3 = static_cast<test_init*>(l2);

    l3->setupEigenfunctionTest(mesh, options, problem, fields);

    /* Now we have a class with testing, which is still really an element :) */
    test_elements.emplace_back(static_cast<test_insert*>(l3));

  }

  PetscReal cycle_time = 1.0; PetscReal max_error = 0;
  RealVec element_error(test_elements.size()); PetscScalar time = 0;
  while (true) {

    std::tie(test_elements, fields) = problem->assembleIntoGlobalDof(std::move(test_elements),
                                                                     std::move(fields),
                                                                     0,
                                                                     mesh->DistributedMesh(),
                                                                     mesh->MeshSection(),
                                                                     options);

    fields = problem->applyInverseMassMatrix(std::move(fields));
    std::tie(fields, time) = problem->takeTimeStep
        (std::move(fields), time, options);

    PetscInt i = 0;
    for (auto &elm: test_elements) {
      auto validate = static_cast<test_init*>(static_cast<test_insert*>(elm.get()));
      element_error(i++) = validate->checkEigenfunctionTestNew(
          mesh, options, time, problem, fields);
    }

    max_error = element_error.maxCoeff() > max_error ? element_error.maxCoeff() : max_error;
    if (time > cycle_time) break;

  }

  PetscReal regression_error = 0.00171; PetscScalar eps = 0.01;
  REQUIRE(max_error <= regression_error * (1 + eps));

}
