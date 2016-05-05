#include <iostream>
#include <Physics/AcousticNew.h>
#include <Element/Element.h>
#include <Element/HyperCube/QuadNew.h>
#include <Element/HyperCube/Quad/QuadP1.h>
#include "catch.h"

TEST_CASE("test_templates", "[template]") {

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
      "--element_shape", "quad",
      "--physics_system", "acoustic",
      "--polynomial_order", "4", NULL};

  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(&argc, &argv, NULL);

  Options options;
  options.setOptions();

  ExodusModel *model = new ExodusModel(options);
  model->initializeParallel();

  Mesh *msh = Mesh::factory(options);
  msh->read(options);
  msh->setupGlobalDof(1, 3, 9, 0, 2, model);

  std::cout << "Hello world." << std::endl;
  auto elm = std::make_shared<QuadNew<AcousticNew<QuadP1>>>(options);

  auto receivers = Receiver::factory(options);
  auto sources = Source::factory(options);

  elm->SetNum(0);
  elm->attachVertexCoordinates(msh->DistributedMesh());
  elm->attachReceiver(receivers);
  elm->attachSource(sources);

  Eigen::MatrixXd dummy(1,1);
  double detJ;
  Eigen::Matrix2d jInv;
  std::tie(jInv, detJ) = QuadP1::inverseJacobianAtPoint(0, 0, elm->VtxCrd());
//  std::cout << elm->computeStiffnessTerm(dummy);




}