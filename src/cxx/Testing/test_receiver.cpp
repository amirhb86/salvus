#include "catch.h"
#include <Element/ElementNew.h>
#include <Element/Element.h>
#include <Problem/Problem.h>

/**
 * TODO
 * This test only runs with mpirun = 1. I think this has something to do with the re-initializaiton
 * of the exodus file. It might be an issue with Petsc never calling ex_close, but i'm not sure.
 */
TEST_CASE("test_receiver", "[receiver]") {
// Set up custom command line arguments.
  PetscOptionsClear();
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--duration", "0.7071067811865475",
      "--time_step", "0.003",
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
      "--receiver_names", "rec1,rec2",
      "--receiver_file_name", "test_rec.h5",
      "--exodus_file_name", "../../salvus_data/unit_test_meshes/homogeneous_iso_cartesian_2D_50s.e",
      "--exodus_model_file_name", "../../salvus_data/unit_test_meshes/homogeneous_iso_cartesian_2D_50s.e",
      "--mesh_type", "newmark",
      "--element_shape", "quad_new",
      "--physics_system", "acoustic",
      "--polynomial_order", "4", NULL};

  int max_order = 4;
  for (int i = 4; i <= max_order; i++) {
    char **argv = const_cast<char **> (arg);
    int argc = sizeof(arg) / sizeof(const char *) - 1;

    PetscOptionsInsert(&argc, &argv, NULL);
    PetscOptionsSetValue("--polynomial_order", std::to_string(i).c_str());

    Options options;
    options.setOptions();

    Mesh *mesh = Mesh::factory(options);
    mesh->read(options);

    ExodusModel *model = new ExodusModel(options);
    model->initializeParallel();

    std::shared_ptr<ElementNew> reference_element = ElementNew::Factory(options);

    Problem *problem = Problem::factory(options.ProblemType());
    problem->initialize(mesh, model, reference_element, options);
    problem->solve(options);

    delete problem;

    if (i > 4) { break; }
  }

}
