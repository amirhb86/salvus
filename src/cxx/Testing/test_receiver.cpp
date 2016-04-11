#include "catch.h"
#include <Element/Element.h>

TEST_CASE("test_receiver", "[receiver]") {
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
      "--number_of_receivers", "2",
      "--receiver_location_x1", "26000,76000",
      "--receiver_location_x2", "26000,76000",
      "--exodus_file_name", "homogeneous_iso_cartesian_2D_50s.e",
      "--exodus_model_file_name", "homogeneous_iso_cartesian_2D_50s.e",
      "--mesh_type", "newmark",
      "--element_shape", "quad",
      "--physics_system", "acoustic",
      "--polynomial_order", "4", NULL};

  int max_order = 10;
  for (int i = 1; i <= max_order; i++) {
    char **argv = const_cast<char **> (arg);
    int argc = sizeof(arg) / sizeof(const char *) - 1;

    PetscOptionsInsert(&argc, &argv, NULL);
    PetscOptionsSetValue("--polynomial_order", std::to_string(i).c_str());

    Options options;
    options.setOptions();

    auto sources = Source::factory(options);
    auto receivers = Receiver::factory(options);
    Mesh *msh = Mesh::factory(options);
    msh->read(options);
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

      // Attach receiver.
      element->attachReceiver(receivers);

    }
  }
}
