#include "catch.h"
#include <Element/Element.h>
#include <Problem/Problem.h>

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
      "--exodus_file_name", "../../salvus_data/unit_test_meshes/homogeneous_iso_cartesian_2D_50s.e",
      "--exodus_model_file_name", "../../salvus_data/unit_test_meshes/homogeneous_iso_cartesian_2D_50s.e",
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

    Mesh *mesh = Mesh::factory(options);
    mesh->read(options);

    ExodusModel *model = new ExodusModel(options);
    model->initializeParallel();

    std::shared_ptr<Element> reference_element = Element::factory(options);

    Problem *problem = Problem::factory(options.ProblemType());
    problem->initialize(mesh, model, reference_element, options);
    problem->solve(options);

    delete problem;

    if (i > 4) { break; }
  }

}
