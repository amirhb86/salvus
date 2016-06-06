#pragma once

// stl.
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
  inline PetscBool SaveMovie() { return mSaveMovie; }
  inline PetscBool TestIC() { return mTestIC; }
  inline PetscBool DisplayDiagnostics() { return mDisplayDiagnostics; }
  
  // Integer getters
  inline PetscInt Dimension() { return mDimension; }
  inline PetscInt PolynomialOrder() { return mPolynomialOrder; }
  inline PetscInt NumberSources() { return mNumberSources; }
  inline PetscInt NumberReceivers() { return mNumRec; }
  inline PetscInt SaveFrameEvery() { return mSaveFrameEvery; }
  inline PetscInt DisplayDiagnosticsEvery() { return mDisplayDiagnosticsEvery; }
  
  
  // Double getters
  inline double Duration() { return mDuration; }
  inline double TimeStep() { return mTimeStep; }
  inline PetscReal IC_Center_x() { return mCenter_x; }
  inline PetscReal IC_Center_y() { return mCenter_y; }
  inline PetscReal IC_Center_z() { return mCenter_z; }
  inline PetscReal IC_SquareSide_L() { return mSquareSide_L; }

  // String getters
  inline std::string PhysicsSystem() { return mPhysicsSystem; }
  inline std::string ExodusMeshFile() { return mExodusMeshFile; }
  inline std::string MeshType() { return mMeshType; }
  inline std::string ElementShape() { return mElementShape; }
  inline std::string ExodusModelFile() { return mExodusModelFile; }
  inline std::string SourceType() { return mSourceType; }
  inline std::string OutputMovieFile() { return mOutputMovieFile; }
  inline std::string ProblemType() { return mProblemType; }
  inline std::string ReceiverType() { return "hdf5"; } // TODO: GENERAL
  inline std::string ReceiverFileName() { return mReceiverFileName; }

  inline std::vector<std::string> DirichletBoundaries() { return mDirichletBoundaryNames; }
  inline std::vector<std::string> RecNames() { return mRecNames; }

  // Vector getters
  inline std::vector<double> RecLocX1() { return mRecLocX1; }
  inline std::vector<double> RecLocX2() { return mRecLocX2; }
  inline std::vector<double> RecLocX3() { return mRecLocX3; }
  inline std::vector<double> SourceLocationX() { return mSourceLocationX; }
  inline std::vector<double> SourceLocationY() { return mSourceLocationY; }
  inline std::vector<double> SourceLocationZ() { return mSourceLocationZ; }
  inline std::vector<double> SourceRickerAmplitude() { return mSourceRickerAmplitude; }
  inline std::vector<double> SourceRickerCenterFreq() { return mSourceRickerCenterFreq; }
  inline std::vector<double> SourceRickerTimeDelay() { return mSourceRickerTimeDelay; }

  // for setting timestep automatically.
  void SetTimeStep(double timestep);
  
  // Setters for testing.
  inline void __SetPolynomialOrder(int ord) { mPolynomialOrder = ord; }
  inline void __SetSourceType(std::string type) { mSourceType = type; }

};
