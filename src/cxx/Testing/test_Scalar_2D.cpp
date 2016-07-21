#include <iostream>
#include <Utilities/Types.h>
#include <Mesh/Mesh.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Utilities/Logging.h>
#include <Physics/Scalar.h>
#include <Element/HyperCube/QuadP1.h>
#include <Element/Simplex/TriP1.h>
#include <Element/ElementAdapter.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Element/Simplex/Triangle.h>
#include <Problem/Problem.h>
#include <Problem/Problem.h>
#include <petscviewerhdf5.h>
#include "catch.h"


template <typename Element>
class TestPlugin: public Element {

 public:
  TestPlugin<Element>(std::unique_ptr<Options> const &options): Element(options) {};

  void setupEigenfunctionTest(std::unique_ptr<Mesh> const &mesh,
                              std::unique_ptr<Options> const &options,
                              std::unique_ptr<Problem> &problem,
                              FieldDict &fields) {

    /* This is hardcoded for the unit test mesh. */
    PetscScalar x0 = 5e4, y0 = 5e4, L = 1e5;
    RealVec pts_x, pts_y;
    std::tie(pts_x, pts_y) = Element::buildNodalPoints();
    RealVec un = (M_PI / L * (pts_x.array() - (x0 + L / 2))).sin() *
        (M_PI / L * (pts_y.array() - (y0 + L / 2))).sin();
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

    RealVec pts_x, pts_y;
    PetscReal x0, y0, L;
    RealVec exact;
    // create `exact` for 2D. We created 2D meshes for
    // quads/triangles with consistent dimensions.
    
    x0 = 5e4; y0 = 5e4; L = 1e5;
    std::tie(pts_x, pts_y) = Element::buildNodalPoints();
    RealVec un_xy = (M_PI / L * (pts_x.array() - (x0 + L / 2))).sin() *
      (M_PI / L * (pts_y.array() - (y0 + L / 2))).sin();
    PetscReal vp = Element::ParAtIntPts("VP").mean();
    PetscReal un_t = cos(M_PI / L * sqrt(2) * time * vp);
    exact = un_t * un_xy;
    
    RealVec u = problem->getFieldOnElement(
        "u", Element::ElmNum(), Element::ClsMap(),
        mesh->DistributedMesh(), mesh->MeshSection(), fields);

    PetscReal element_error = (exact - u).array().abs().maxCoeff();
    return element_error;

  }
};

typedef ElementAdapter<TestPlugin<Scalar<TensorQuad<QuadP1>>>> test_insert_quadP1;
typedef TestPlugin<Scalar<TensorQuad<QuadP1>>> test_init_quadP1;
typedef ElementAdapter<Scalar<TensorQuad<QuadP1>>> unguard_quadP1;
typedef Scalar<TensorQuad<QuadP1>> raw_quadP1;

typedef ElementAdapter<TestPlugin<Scalar<Triangle<TriP1>>>> test_insert_triP1;
typedef TestPlugin<Scalar<Triangle<TriP1>>> test_init_triP1;
typedef ElementAdapter<Scalar<Triangle<TriP1>>> unguard_triP1;
typedef Scalar<Triangle<TriP1>> raw_triP1;

PetscReal runEigenFunctionTestQuad(std::vector<std::unique_ptr<Element>> test_elements,
                               std::unique_ptr<Mesh> &mesh,
                               std::unique_ptr<ExodusModel> const &model,
                               std::unique_ptr<Options> const &options,
                               std::unique_ptr<Problem> &problem,
                               FieldDict &fields,
                               PetscReal cycle_time
                               ) {

  RealVec element_error(test_elements.size()); PetscScalar time = 0;
  PetscReal max_error = 0.0;
  while (true) {

    std::tie(test_elements, fields) = problem->assembleIntoGlobalDof(
        std::move(test_elements), std::move(fields),
        mesh->DistributedMesh(), mesh->MeshSection(),
        options);

    fields = problem->applyInverseMassMatrix(std::move(fields));
    std::tie(fields, time) = problem->takeTimeStep
        (std::move(fields), time, options);

    PetscInt i = 0;
    for (auto &elm: test_elements) {            
        auto validate = static_cast<test_init_quadP1*>(static_cast<test_insert_quadP1*>(elm.get()));
        element_error(i++) = validate->checkEigenfunctionTestNew(mesh, options, time, problem, fields);
    }

    // only saves if '--saveMovie' is set in command line options
    problem->saveSolution(time, {"u"}, fields, mesh->DistributedMesh());

    std::cout << "TIME:      " << time << std::endl;
    max_error = element_error.maxCoeff() > max_error ? element_error.maxCoeff() : max_error;
    if (time > cycle_time) break;
    
  }
  return max_error;  
}

PetscReal runEigenFunctionTest(std::vector<std::unique_ptr<Element>> test_elements,
                               std::unique_ptr<Mesh> &mesh,
                               std::unique_ptr<ExodusModel> const &model,
                               std::unique_ptr<Options> const &options,
                               std::unique_ptr<Problem> &problem,
                               FieldDict &fields,
                               PetscReal cycle_time,
                               ElementType element_type
                               ) {

  RealVec element_error(test_elements.size()); PetscScalar time = 0;
  PetscReal max_error = 0.0;
  while (true) {

    std::tie(test_elements, fields) = problem->assembleIntoGlobalDof(
        std::move(test_elements), std::move(fields),
        mesh->DistributedMesh(), mesh->MeshSection(),
        options);

    fields = problem->applyInverseMassMatrix(std::move(fields));
    std::tie(fields, time) = problem->takeTimeStep
        (std::move(fields), time, options);

    PetscInt i = 0;
    for (auto &elm: test_elements) {
      switch(element_type) {
      case ElementType::QUADP1:
        {
          auto validate = static_cast<test_init_quadP1*>(static_cast<test_insert_quadP1*>(elm.get()));
          element_error(i++) = validate->checkEigenfunctionTestNew(mesh, options, time, problem, fields);
          break;
        }
      case ElementType::TRIP1:
        {
          auto validate = static_cast<test_init_triP1*>(static_cast<test_insert_triP1*>(elm.get()));
          element_error(i++) = validate->checkEigenfunctionTestNew(mesh, options, time, problem, fields);
          break;
        }
      default:
        ERROR() << "Element type not supported";
      }
    }

    // only saves if '--saveMovie' is set in command line options
    problem->saveSolution(time, {"u"}, fields, mesh->DistributedMesh());

    std::cout << "TIME:      " << time << std::endl;
    max_error = element_error.maxCoeff() > max_error ? element_error.maxCoeff() : max_error;
    if (time > cycle_time) break;
    
  }
  return max_error;  
}


TEST_CASE("Test analytic eigenfunction solution for scalar "
              "equation in 2D with quadrilateral", "[quad_eigenfunction]") {

  std::string e_file = "quad_eigenfunction.e";

  PetscOptionsClear(NULL);
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--element_shape", "quad_new",
      "--exodus_file_name", e_file.c_str(),
      "--exodus_model_file_name", e_file.c_str(),
      "--polynomial_order", "4", NULL};
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(NULL, &argc, &argv, NULL);

  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  std::unique_ptr<Problem> problem(Problem::Factory(options));
  std::unique_ptr<ExodusModel> model(new ExodusModel(options));
  std::unique_ptr<Mesh> mesh(Mesh::Factory(options));

  model->initializeParallel();
  mesh->read(options);
  mesh->setupGlobalDof(1, 3, 9, 0, 2, model);
  
  PetscReal cycle_time = 24.39; PetscReal max_error = 0;
  cycle_time = 0.05;

  std::vector<std::unique_ptr<Element>> test_elements;
  auto elements = problem->initializeElements(mesh, model, options);
  auto fields = problem->initializeGlobalDofs(elements, mesh);

  /* Rip apart elements and insert testing mixin. */
  for (auto &e: elements) {

    /* Rip out the master Element class. */
    auto l1 = static_cast<unguard_quadP1*>(e.release());

    /* Rip out the Element adapter. */
    auto l2 = static_cast<raw_quadP1*>(l1);

    /* Attach the tester. */
    auto l3 = static_cast<test_init_quadP1*>(l2);

    l3->setupEigenfunctionTest(mesh, options, problem, fields);

    /* Now we have a class with testing, which is still really an element :) */
    test_elements.emplace_back(static_cast<test_insert_quadP1*>(l3));

  }

  runEigenFunctionTest(std::move(test_elements),mesh,model,options,problem,fields,cycle_time, ElementType::QUADP1);
  
  PetscReal regression_error = 0.001288; PetscScalar eps = 0.01;
  REQUIRE(max_error <= regression_error * (1 + eps));

}

