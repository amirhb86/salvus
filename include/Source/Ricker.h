#pragma once

// stl.
#include <vector>
#include <memory>

// 3rd party.
#include <petsc.h>

// parents.
#include <Source/Source.h>

// forward decl.
class Options;

class Ricker: public Source {

  double mAmplitude;
  double mTimeDelay;
  double mCenterFreq;

  Eigen::VectorXd mDirection;

 public:

  Ricker(std::unique_ptr<Options> const &options);
  ~Ricker() {};
  Eigen::VectorXd fire(const double &time, const PetscInt &time_idx);

};

