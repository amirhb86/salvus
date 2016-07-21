#include <petsc.h>
#include <Utilities/Options.h>
#include <Source/Source.h>
#include <Mesh/Mesh.h>
#include <Model/ExodusModel.h>
#include <Problem/Problem.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>
#include <Physics/Scalar.h>
#include "catch.h"

using namespace std;

double true_ricker(const PetscReal time, const PetscReal freq,
                   const PetscReal delay, const PetscReal amp) {
  double factor = M_PI * M_PI * freq * freq * (time - delay) * (time - delay);
  return amp * ((1 - 2 * factor) * exp(-1 * factor));
}

/* Test plugins. */
template <typename Element>
class QuadTestPlugin: public Element {
 public:
  QuadTestPlugin<Element>(unique_ptr<Options> const &options): Element(options) {};
};

typedef class Scalar<TensorQuad<QuadP1>> quadtest_p1;


TEST_CASE("Test source functionality", "[source]") {

  SECTION("unit") {

    PetscOptionsClear(NULL);
    const char *arg[] =
        {"salvus_test", "--testing", "true", "--number_of_sources", "2", "--source_type", "ricker",
         "--source_location_x", "50000,50000", "--source_location_y", "0,0", "--source_location_z",
         "80000,90000", "--ricker_amplitude", "10,20", "--ricker_time_delay", "0.1,0.01",
         "--ricker_center_freq", "50,60", NULL};

    char **argv = const_cast<char **> (arg);
    int argc = sizeof(arg) / sizeof(const char *) - 1;
    PetscOptionsInsert(NULL, &argc, &argv, NULL);

    SECTION("general") {

      vector<PetscReal> x{50000, 50000};
      vector<PetscReal> y{0, 0};
      vector<PetscReal> z{80000, 90000};

      /* Need something to complete the source, so we choose ricker. */
      PetscOptionsSetValue(NULL, "--ricker_amplitude", "10,20");
      PetscOptionsSetValue(NULL, "--ricker_time_delay", "0.1,0.01");
      PetscOptionsSetValue(NULL, "--ricker_center_freq", "50,60");
      std::unique_ptr<Options> options(new Options);
      options->setOptions();
      auto sources = Source::Factory(options);

      /* Require that naming worked correctly. */
      REQUIRE(Source::NumSources() == 2);

      for (PetscInt i = 0; i < Source::NumSources(); i++) {
        REQUIRE(sources[i]->LocX() == x[i]);
        REQUIRE(sources[i]->LocY() == y[i]);
        REQUIRE(sources[i]->LocZ() == z[i]);
      }

      /* Require deletion works properly. */
      sources.back().reset();
      REQUIRE(Source::NumSources() == 1);
      sources[0].reset();
      REQUIRE(Source::NumSources() == 0);

    }

    SECTION("ricker") {

      /* Options for ricker. */
      PetscOptionsSetValue(NULL, "--ricker_amplitude", "10,20");
      PetscOptionsSetValue(NULL, "--ricker_time_delay", "0.1,0.01");
      PetscOptionsSetValue(NULL, "--ricker_center_freq", "50,60");
      std::unique_ptr<Options> options(new Options);
      options->setOptions();

      /* True options. */
      vector<PetscReal> ricker_amp{10, 20};
      vector<PetscReal> ricker_time{0.1, 0.01};
      vector<PetscReal> ricker_freq{50, 60};

      auto sources = Source::Factory(options);

      /* Require ricker source fires properly. */
      for (PetscInt i = 0; i < Source::NumSources(); i++) {
        REQUIRE(sources[i]->Num() == i);
        for (PetscInt j = 0; j < 1000; j++) {
          PetscReal time = j * 1e-3;
          REQUIRE(sources[i]->fire(time)
                      == true_ricker(time, ricker_freq[i], ricker_time[i], ricker_amp[i]));
        }
      }

    }

    SECTION("integration") {

      string e_file = "fluid_layer_over_elastic_cartesian_2D_50s.e";
      PetscOptionsSetValue(NULL, "--exodus_file_name", e_file.c_str());
      PetscOptionsSetValue(NULL, "--exodus_model_file_name", e_file.c_str());
      PetscOptionsSetValue(NULL, "--element_shape", "quad_new");
      PetscOptionsSetValue(NULL, "--polynomial_order", "9");
      std::unique_ptr<Options> options(new Options);
      options->setOptions();

      unique_ptr<ExodusModel> model(new ExodusModel(options));
      model->initializeParallel();

      unique_ptr<Mesh> mesh(new Mesh(options));
      mesh->read(options);
      mesh->setupGlobalDof(2, model, options);

      /* True values. */
      vector<int> src_elm{1, 0, 1, 0};
      vector<PetscReal> ricker_amp{10, 20};
      vector<PetscReal> ricker_time{0.1, 0.01};

      auto problem = Problem::Factory(options);
      auto elements = problem->initializeElements(mesh, model, options);

      SECTION("quad") {
        for (auto &e: elements) {
          if (e->Num() == 0) {
            REQUIRE(e->computeSourceTerm(ricker_time[0]).sum() == Approx(ricker_amp[0]));
          }
          else if (e->Num() == 2) {
            REQUIRE(e->computeSourceTerm(ricker_time[1]).sum() == Approx(ricker_amp[1]));
          }
          else {
            REQUIRE(e->computeSourceTerm(ricker_time[0]).sum() == Approx(0.0));
            REQUIRE(e->computeSourceTerm(ricker_time[1]).sum() == Approx(0.0));
          }
        }
      }
    }
  }
}
