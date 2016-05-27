#include "catch.h"
#include <Element/Element.h>
#include <Element/ElementAdapter.h>
#include <Element/HyperCube/Quad.h>
#include <Element/HyperCube/QuadP1.h>
#include <Source/Source.h>
#include <Mesh/Mesh.h>
#include <Model/ExodusModel.h>

TEST_CASE("test quadbase", "[quad]") {


  SECTION("Sources") {

    int order = 4;

    // Set up custom command line arguments.
    PetscOptionsClear();
    const char *arg[] = {
        "salvus_test",
        "--testing", "true",
        "--number_of_sources", "1",
        "--source_type", "ricker",
        "--source_location_x", "76000",
        "--source_location_z", "76000",
        "--ricker_amplitude", "10.0",
        "--ricker_time_delay", "1.0",
        "--ricker_center_freq", "1.0",
        "--exodus_file_name", "homogeneous_iso_cartesian_2D_50s.e",
        "--exodus_model_file_name", "homogeneous_iso_cartesian_2D_50s.e",
        "--mesh_type", "newmark",
        "--element_shape", "quad_new",
        "--physics_system", "acoustic",
        "--polynomial_order", "4", NULL};

    char **argv = const_cast<char **> (arg);
    int argc = sizeof(arg) / sizeof(const char *) - 1;

    PetscOptionsInsert(&argc, &argv, NULL);
    PetscOptionsSetValue("--polynomial_order", std::to_string(order).c_str());

    Options options;
    options.setOptions();

    auto sources = Source::factory(options);
    Mesh *msh = Mesh::factory(options);
    msh->read(options);
    std::shared_ptr<Element> elm = Element::Factory({"u"}, {}, options);
    std::vector<std::shared_ptr<Element>> elms;

    ExodusModel *model = new ExodusModel(options);
    model->initializeParallel();

    msh->setupGlobalDof(elm->NumDofVtx(), elm->NumDofEdg(), elm->NumDofFac(),
                        elm->NumDofVol(), elm->NumDim(), model);

    // Get a list of all local elements.
    for (int i = 0; i < msh->NumberElementsLocal(); i++) {
      elms.push_back(elm->clone());
    }

    // Set up elements.
    int element_number = 0;
    for (auto &element : elms) {

      // Give each element a number starting from zero.
      element->SetNum(element_number++);

      // Get vertex coordinates from the PETSc DMPLEX.
      element->attachVertexCoordinates(msh->DistributedMesh());

      // Attach source.
      element->attachSource(sources);

    }

    SECTION("Source interpolation is correct.") {

      Eigen::VectorXd regression_source(25);
      regression_source <<
          3.29303e-10, -2.11293e-10, 3.21071e-09, 2.38793e-10, -3.56745e-10,
          -2.11293e-10, 1.35573e-10, -2.0601e-09, -1.53218e-10, 2.289e-10,
          3.21071e-09, -2.0601e-09, 3.13044e-08, 2.32824e-09, -3.47827e-09,
          2.38793e-10, -1.53218e-10, 2.32824e-09, 1.7316e-10, -2.58693e-10,
          -3.56745e-10, 2.289e-10, -3.47827e-09, -2.58693e-10, 3.86474e-10;

      for (auto &element : elms) {
        Eigen::MatrixXd test = element->computeSourceTerm(1.0);
        if (test.maxCoeff() > 0) {
          REQUIRE(test.isApprox(regression_source, 1e-6));
        }
      }
    }

  }

}