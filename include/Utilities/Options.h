#pragma once

// stl.
#include <iostream>
#include <iosfwd>
#include <string>
#include <vector>

// 3rd party.
#include <petsc.h>

class Options {

  // Integer options.
  PetscInt mDimension;
  PetscInt mNumberSources;
  PetscInt mPolynomialOrder;

  // Double options.
  // Simulation duration
  double mDuration;
  // Time step.
  double mTimeStep;

  // String options.
  std::string mMeshType;
  std::string mExodusMeshFile;
  std::string mElementShape;
  std::string mPhysicsSystem;
  std::string mExodusModelFile;
  std::string mSourceType;
  std::string mProblemType;
  std::string mTimeStepType;

  std::vector<std::string> mDirichletBoundaryNames;
  std::vector<std::string> mRecNames;

  // Determines the output name of the movie.
  std::string mOutputMovieFile;
  // Save movie?, and how often.
  PetscBool mSaveMovie;
  PetscInt mSaveFrameEvery;

  // Diagnostic information
  PetscBool mDisplayDiagnostics;
  PetscInt mDisplayDiagnosticsEvery;
  
  // Run initial condition test with exact solution
  PetscBool mTestIC;
  // exact solution parameters
  PetscReal mCenter_x;
  PetscReal mCenter_y;
  PetscReal mCenter_z;
  PetscReal mSquareSide_L;

  // Sources.
  // TODO: Move to simple HDF5 file.
  std::vector<double> mSourceLocationX;
  std::vector<double> mSourceLocationY;
  std::vector<double> mSourceLocationZ;
  std::vector<double> mSourceRickerAmplitude;
  std::vector<double> mSourceRickerCenterFreq;
  std::vector<double> mSourceRickerTimeDelay;

  // Receivers.
  // TODO: Move to simple HDF5 file.
  PetscInt mNumRec;
  std::string mReceiverFileName;
  std::vector<double> mRecLocX1;
  std::vector<double> mRecLocX2;
  std::vector<double> mRecLocX3;

 public:

  PetscErrorCode setOptions();

  // Bool getters
  inline PetscBool SaveMovie() const { return mSaveMovie; }
  inline PetscBool TestIC() const { return mTestIC; }
  inline PetscBool DisplayDiagnostics() const { return mDisplayDiagnostics; }
  
  // Integer getters
  inline PetscInt Dimension() const { return mDimension; }
  inline PetscInt PolynomialOrder() const { return mPolynomialOrder; }
  inline PetscInt NumberSources() {  return mNumberSources; }
  inline PetscInt NumberReceivers() const { return mNumRec; }
  inline PetscInt SaveFrameEvery() const { return mSaveFrameEvery; }
  inline PetscInt DisplayDiagnosticsEvery() const { return mDisplayDiagnosticsEvery; }
  
  
  // Double getters
  inline double Duration() const { return mDuration; }
  inline double TimeStep() const { return mTimeStep; }
  inline PetscReal IC_Center_x() const { return mCenter_x; }
  inline PetscReal IC_Center_y() const { return mCenter_y; }
  inline PetscReal IC_Center_z() const { return mCenter_z; }
  inline PetscReal IC_SquareSide_L() const { return mSquareSide_L; }

  // String getters
  inline std::string PhysicsSystem() const { return mPhysicsSystem; }
  inline std::string ExodusMeshFile() const { return mExodusMeshFile; }
  inline std::string MeshType() const { return mMeshType; }
  inline std::string ElementShape() const { return mElementShape; }
  inline std::string ExodusModelFile() const { return mExodusModelFile; }
  inline std::string SourceType() const { return mSourceType; }
  inline std::string OutputMovieFile() const { return mOutputMovieFile; }
  inline std::string ProblemType() const { return mProblemType; }
  inline std::string ReceiverType() const { return "hdf5"; } // TODO: GENERAL
  inline std::string ReceiverFileName() const { return mReceiverFileName; }

  inline std::vector<std::string> DirichletBoundaries() const { return mDirichletBoundaryNames; }
  inline std::vector<std::string> RecNames() const { return mRecNames; }

  // Vector getters
  inline std::vector<double> RecLocX1() const { return mRecLocX1; }
  inline std::vector<double> RecLocX2() const { return mRecLocX2; }
  inline std::vector<double> RecLocX3() const { return mRecLocX3; }
  inline std::vector<double> SourceLocationX() const { return mSourceLocationX; }
  inline std::vector<double> SourceLocationY() const { return mSourceLocationY; }
  inline std::vector<double> SourceLocationZ() const { return mSourceLocationZ; }
  inline std::vector<double> SourceRickerAmplitude() const { return mSourceRickerAmplitude; }
  inline std::vector<double> SourceRickerCenterFreq() const { return mSourceRickerCenterFreq; }
  inline std::vector<double> SourceRickerTimeDelay() const { return mSourceRickerTimeDelay; }

  // for setting timestep automatically.
  void SetTimeStep(double timestep);
  
  // Setters for testing.
  inline void __SetPolynomialOrder(int ord) { mPolynomialOrder = ord; }
  inline void __SetSourceType(std::string type) { mSourceType = type; }

};
