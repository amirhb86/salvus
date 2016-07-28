#include <mpi.h>
#include <Source/Source.h>
#include <Utilities/Options.h>
#include <stdexcept>
#include <iostream>
#include <Utilities/Logging.h>

/* Initialize counter. */
PetscInt Source::number = 0;

std::vector<std::unique_ptr<Source>> Source::Factory(std::unique_ptr<Options> const &options) {

  std::vector<std::unique_ptr<Source>> sources;
  for (PetscInt i = 0; i < options->NumberSources(); i++) {
    try {
      if (options->SourceType() == "ricker") {
        sources.push_back(std::unique_ptr<Ricker>(new Ricker(options)));
      } else {
        throw std::runtime_error("Runtime error: Source type " + options->SourceType() + " not supported.");
      }

    } catch (std::exception &e) {
      LOG() << e.what();
      MPI_Abort(PETSC_COMM_WORLD, -1);
    }
  }

  return sources;
}

Source::Source(std::unique_ptr<Options> const &options) {

  /* Increment global number, save this particular one. */
  SetNum(number++);

  /* Set locations. */
  SetLocX(options->SrcLocX()[Num()]);
  SetLocY(options->SrcLocY()[Num()]);
  if (options->SrcLocZ().size()) {
    SetLocZ(options->SrcLocZ()[Num()]);
  } else {
    SetLocZ(0.0);
  }

}

Source::~Source() { --number; }

Ricker::Ricker(std::unique_ptr<Options> const &options): Source(options) {

  mTimeDelay = options->SrcRickerTimeDelay()[Num()];
  mAmplitude = options->SrcRickerAmplitude()[Num()];
  mCenterFreq = options->SrcRickerCenterFreq()[Num()];

}

double Ricker::fire(const double &time) {

  double factor = M_PI * M_PI * mCenterFreq * mCenterFreq * (time - mTimeDelay) * (time - mTimeDelay);
  return mAmplitude * ((1 - 2 * factor) * exp(-1 * factor));

}
