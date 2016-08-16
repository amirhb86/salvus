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

 public:

  Ricker(std::unique_ptr<Options> const &options);
  ~Ricker() {};
  double fire(const double &time);

};

