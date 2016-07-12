#include <iostream>
#include <Utilities/Types.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Physics/Scalar.h>
#include <Element/HyperCube/QuadP1.h>
#include <Element/ElementAdapter.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Problem/ProblemNew.h>
#include <Mesh/ElasticAcousticNewmark3D.h>
#include <Problem/ProblemNew.h>
#include <petscviewerhdf5.h>
#include "catch.h"

template <typename Element>
class TestPlugin: public Element {

 public:
  TestPlugin<Element>(std::unique_ptr<Options> const &options): Element(options) {};

  void setupEigenfunctionTest(std::unique_ptr<Mesh> const &mesh,
                              std::unique_ptr<Options> const &options,
                              std::unique_ptr<ProblemNew> &problem,
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

  }
};

typedef ElementAdapter<TestPlugin<Scalar<TensorQuad<QuadP1>>>> test_insert;
typedef TestPlugin<Scalar<TensorQuad<QuadP1>>> test_init;
typedef ElementAdapter<Scalar<TensorQuad<QuadP1>>> unguard;
typedef Scalar<TensorQuad<QuadP1>> raw;

TEST_CASE("Test analytic eigenfunction solution for scalar "
              "equation in 2D", "[quad_eigenfunction]") {

  std::string e_file = "quad_eigenfunction.e";

  PetscOptionsClear();
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--element_shape", "quad_new",
      "--exodus_file_name", e_file.c_str(),
      "--exodus_model_file_name", e_file.c_str(),
      "--polynomial_order", "4", NULL};
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(&argc, &argv, NULL);

  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  std::unique_ptr<ProblemNew> problem(ProblemNew::Factory(options));
  std::unique_ptr<ExodusModel> model(new ExodusModel(options));
  std::unique_ptr<Mesh> mesh(new ElasticAcousticNewmark3D(options));

  model->initializeParallel();
  mesh->read(options);
  mesh->setupGlobalDof(1, 3, 9, 0, 2, model);

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



//  std::string filename = "test.h5";
//  PetscViewer viewer = nullptr;
//  PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(),
//                      FILE_MODE_WRITE, &viewer);
//  PetscViewerHDF5PushGroup(viewer, "/");
//  DMView(mesh->DistributedMesh(), viewer);
//  VecView(fields["u"]->mGlb, viewer);
//  PetscViewerDestroy(&viewer);

  while (true) {

    std::tie(test_elements, fields) = problem->assembleIntoGlobalDof(
        std::move(test_elements), std::move(fields),
        mesh->DistributedMesh(), mesh->MeshSection(),
        options);

    fields = problem->applyInverseMassMatrix(std::move(fields));
    fields = problem->takeTimeStep(std::move(fields));

  }

}