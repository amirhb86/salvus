#include <Utilities/Options.h>
#include <Source/Source.h>
#include <Mesh/Mesh.h>
#include <Mesh/ElasticAcousticNewmark3D.h>
#include <Model/ExodusModel.h>
#include <Problem/ProblemNew.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>
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

typedef class TensorQuad<QuadP1> quadtest_p1;


TEST_CASE("Test source functionality", "[source]") {

  SECTION("unit") {

    PetscOptionsClear();
    const char *arg[] = {
        "salvus_test",
        "--testing", "true",
        "--number_of_sources", "2",
        "--source_type", "ricker",
        "--source_location_x", "50000,50000",
        "--source_location_y", "0,0",
        "--source_location_z", "80000,90000",
        "--ricker_amplitude", "10,20",
        "--ricker_time_delay", "0.1,0.01",
        "--ricker_center_freq", "50,60",
        NULL
    };

    char **argv = const_cast<char **> (arg);
    int argc = sizeof(arg) / sizeof(const char *) - 1;
    PetscOptionsInsert(&argc, &argv, NULL);

    SECTION("general") {

      vector<PetscReal> x{50000, 50000};
      vector<PetscReal> y{0, 0};
      vector<PetscReal> z{80000, 90000};

      /* Need something to complete the source, so we choose ricker. */
      PetscOptionsSetValue("--ricker_amplitude", "10,20");
      PetscOptionsSetValue("--ricker_time_delay", "0.1,0.01");
      PetscOptionsSetValue("--ricker_center_freq", "50,60");
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
      sources.back().reset(); REQUIRE(Source::NumSources() == 1);
      sources[0].reset(); REQUIRE(Source::NumSources() == 0);

    }

    SECTION("ricker") {

      /* Options for ricker. */
      PetscOptionsSetValue("--ricker_amplitude", "10,20");
      PetscOptionsSetValue("--ricker_time_delay", "0.1,0.01");
      PetscOptionsSetValue("--ricker_center_freq", "50,60");
      std::unique_ptr<Options> options(new Options);
      options->setOptions();

      /* True options. */
      vector<PetscReal> ricker_amp  {10, 20};
      vector<PetscReal> ricker_time {0.1, 0.01};
      vector<PetscReal> ricker_freq {50, 60};

      auto sources = Source::Factory(options);

      /* Require ricker source fires properly. */
      for (PetscInt i = 0; i < Source::NumSources(); i++) {
        REQUIRE(sources[i]->Num() == i);
        for (PetscInt j = 0; j < 1000; j++) {
          PetscReal time = j * 1e-3;
          REQUIRE(sources[i]->fire(time) == true_ricker(time, ricker_freq[i],
                                                        ricker_time[i], ricker_amp[i]));
        }
      }

    }

    SECTION("integration") {

      string e_file = "fluid_layer_over_elastic_cartesian_2D_50s.e";
      PetscOptionsSetValue("--exodus_file_name", e_file.c_str());
      PetscOptionsSetValue("--exodus_model_file_name", e_file.c_str());
      PetscOptionsSetValue("--element_shape", "quad_new");
      PetscOptionsSetValue("--polynomial_order", "1");
      std::unique_ptr<Options> options(new Options);
      options->setOptions();

      unique_ptr<ExodusModel> model(new ExodusModel(options));
      model->initializeParallel();

      unique_ptr<Mesh> mesh(new ElasticAcousticNewmark3D(options));
      mesh->read(options);
      mesh->setupGlobalDof(1, 3, 9, 0, 2, model);

      /* True values. */
      vector<int>       src_elm     {1, 0, 1, 0};
      vector<PetscReal> ricker_amp  {10, 20};
      vector<PetscReal> ricker_time {0.1, 0.01};

      auto problem = ProblemNew::Factory(options);
      auto elements = problem->initializeElements(mesh, model, options);

      SECTION("quad") {
        vector<quadtest_p1*> p1quads;
        for (auto &e: elements) { p1quads.push_back(dynamic_cast<quadtest_p1*> (e.get())); }
        int i = 0;
        for (auto &e: p1quads) {
          /* Ensure sources are on the correct elements. */
          REQUIRE(e->Sources().size() == src_elm[i]);
          for (auto &s: e->Sources()) {
            if (i == 0) {
              REQUIRE(e->integrateField(s->fire(0.1) *
                  e->getDeltaFunctionCoefficients(s->LocR(), s->LocS())) == ricker_amp[0]);
            }
            if (i == 2) {
              REQUIRE(e->integrateField(s->fire(0.01) *
                  e->getDeltaFunctionCoefficients(s->LocR(), s->LocS())) == ricker_amp[1]);
            }
          }
          i++;
        }
      }

    }

  }

}
