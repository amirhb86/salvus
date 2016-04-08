#include "catch.h"
#include <Element/Element.h>
#include <Element/HyperCube/Quad/Acoustic.h>

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
      "--ricker_amplitude", "10.0",
      "--ricker_time_delay", "1.0",
      "--ricker_center_freq", "1.0",
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

    }

    SECTION("Interpolation is correct.") {
      for (auto source : sources) {
        REQUIRE(source->PhysicalLocationX() == Approx(76000));
        REQUIRE(source->PhysicalLocationZ() == Approx(76000));
        REQUIRE(source->ReferenceLocationEps() == Approx(0.04));
        REQUIRE(source->ReferenceLocationEta() == Approx(0.04));
      }
    }

    SECTION("Point source integrates properly.") {
      for (auto &elem : elms) {
        Quad *quad = dynamic_cast<Quad *> (elem);
        if (quad->NumSrc()) {
          REQUIRE(quad->integrateField(quad->computeSourceTerm(1.0)) == Approx(10));
        } else {
          REQUIRE(quad->computeSourceTerm(1.0).isZero());
        }
      }
    }
  }
}