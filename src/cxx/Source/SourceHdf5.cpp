#include <mpi.h>
#include <Source/Source.h>
#include <Source/SourceHdf5.h>
#include <Utilities/Options.h>
#include <stdexcept>
#include <iostream>
#include <Utilities/Logging.h>

SourceHdf5::SourceHdf5(std::unique_ptr<Options> const &options): Source(options) {

  /* Set locations. */
  SetLocX(options->SrcLocX()[Num()]);
  SetLocY(options->SrcLocY()[Num()]);
  if (options->SrcLocZ().size()) {
    SetLocZ(options->SrcLocZ()[Num()]);
  }

  mTimeDelay = options->SrcRickerTimeDelay()[Num()];
  mAmplitude = options->SrcRickerAmplitude()[Num()];
  mCenterFreq = options->SrcRickerCenterFreq()[Num()];

}

double SourceHdf5::fire(const double &time) {

  double factor = M_PI * M_PI * mCenterFreq * mCenterFreq * (time - mTimeDelay) * (time - mTimeDelay);
  return mAmplitude * ((1 - 2 * factor) * exp(-1 * factor));

}
