#include "catch.h"
#include <Element/Element.h>
#include <Element/HyperCube/Quad/Acoustic.h>

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

TEST_CASE("point_source_interpolation", "[source]") {

  // Set up custom command line arguments.
  PetscOptionsClear();
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--number_of_sources", "1",
      "--source_type", "ricker",
      "--source_location_x", "76000",
      "--source_location_z", "76000",
      "--ricker_amplitude", "1.0",
      "--ricker_time_delay", "1.0",
      "--ricker_center_freq", "1.0",
      "--exodus_file_name", "homogeneous_iso_cartesian_2D_50s.e",
      "--exodus_model_file_name", "homogeneous_iso_cartesian_2D_50s.e",
      "--mesh_type", "newmark",
      "--element_shape", "quad",
      "--physics_system", "acoustic",
      "--polynomial_order", "4", NULL};

  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(&argc, &argv, NULL);

  Options options;
  options.setOptions();

  auto sources = Source::factory(options);
  Mesh *msh = Mesh::factory(options); msh->read(options);
  Element *elm = Element::factory(options);
  std::vector<Element *> elms;

  msh->setupGlobalDof(elm->NumDofVtx(), elm->NumDofEdg(),
                      elm->NumDofFac(), elm->NumDofVol(),
                      elm->NumDim());

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

  for (auto source : sources) {
    REQUIRE(source->ReferenceLocationEps() == Approx(0.04));
    REQUIRE(source->ReferenceLocationEta() == Approx(0.04));
  }

}