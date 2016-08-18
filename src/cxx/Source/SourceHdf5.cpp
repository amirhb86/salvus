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

Eigen::VectorXd SourceHdf5::fire(const double &time, const PetscInt &time_idx) {

  Eigen::VectorXd src(mNumComponents);
  return src;

}
