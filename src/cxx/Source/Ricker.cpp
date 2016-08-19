#include <mpi.h>
#include <Source/Source.h>
#include <Source/Ricker.h>
#include <Utilities/Options.h>
#include <stdexcept>
#include <iostream>
#include <Utilities/Logging.h>

Ricker::Ricker(std::unique_ptr<Options> const &options): Source(options) {

  /* Set locations. */
  SetLocX(options->SrcLocX()[Num()]);
  SetLocY(options->SrcLocY()[Num()]);
  if (options->SrcLocZ().size()) {
    SetLocZ(options->SrcLocZ()[Num()]);
  }

  mNumComponents = options->SrcNumComponents()[Num()];

  mTimeDelay = options->SrcRickerTimeDelay()[Num()];
  mAmplitude = options->SrcRickerAmplitude()[Num()];
  mCenterFreq = options->SrcRickerCenterFreq()[Num()];

  mDirection.resize(mNumComponents);
  for (auto i=0; i<mNumComponents; i++)
    mDirection(i) = 0.0;
    mDirection(0) = 1.0;

}

Eigen::VectorXd Ricker::fire(const double &time, const PetscInt &time_idx) {

  const double factor = M_PI * M_PI * mCenterFreq * mCenterFreq * (time - mTimeDelay) * (time - mTimeDelay);
  const double ricker_force = (mAmplitude * ((1 - 2 * factor) * exp(-1 * factor)));
  return ricker_force * mDirection;

}
