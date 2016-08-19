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

class SourceHdf5: public Source {

  double mAmplitude;
  double mTimeDelay;
  double mCenterFreq;

  std::string mSourceFileName;
  std::string mSourceName;

  PetscInt mNumTimeSteps;

 public:

  SourceHdf5(std::unique_ptr<Options> const &options);
  ~SourceHdf5() {};
  Eigen::VectorXd fire(const double &time, const PetscInt &time_idx);
  
  void loadData();

  Eigen::MatrixXd source_time_function;
};