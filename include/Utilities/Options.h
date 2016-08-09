#pragma once

// stl.
#include <iostream>
#include <iosfwd>
#include <string>
#include <vector>

// 3rd party.
#include <petsc.h>

class Options {

  PetscBool mTesting;
  PetscBool mSaveMovie;

  PetscInt mNumDim;
  PetscInt mNumSrc;
  PetscInt mPolynomialOrder;
  PetscInt mSaveFrameEvery;

  PetscReal mDuration;
  PetscReal mTimeStep;

  std::string mMeshFile;
  std::string mModelFile;
  std::string mSourceType;
  std::string mMovieFile;

  // Sources.
  std::vector<PetscReal> mSrcLocX;
  std::vector<PetscReal> mSrcLocY;
  std::vector<PetscReal> mSrcLocZ;
  std::vector<PetscReal> mSrcRickerAmplitude;
  std::vector<PetscReal> mSrcRickerCenterFreq;
  std::vector<PetscReal> mSrcRickerTimeDelay;

  // Receivers.
  PetscInt mNumRec;
  std::string mReceiverFileName;
  std::vector<PetscReal> mRecLocX;
  std::vector<PetscReal> mRecLocY;
  std::vector<PetscReal> mRecLocZ;
  std::vector<std::string> mRecNames;
  std::vector<std::string> mMovieFields;

  // Boundaries.
  std::vector<std::string> mHomogeneousDirichletBoundaries;

 public:

  void setOptions();

  PetscBool SaveMovie() const { return mSaveMovie; }

  PetscInt Dimension() const { return mNumDim; }
  PetscInt PolynomialOrder() const { return mPolynomialOrder; }
  PetscInt NumberSources() {  return mNumSrc; }
  PetscInt NumberReceivers() const { return mNumRec; }

  PetscReal Duration() const { return mDuration; }
  PetscReal TimeStep() const { return mTimeStep; }

  std::string MeshFile() const { return mMeshFile; }
  std::string ReceiverType() const { return "hdf5"; }
  std::string ModelFile() const { return mModelFile; }
  std::string MovieFile() const { return mMovieFile; }
  std::string SourceType() const { return mSourceType; }
  std::string ReceiverFileName() const { return mReceiverFileName; }

  std::vector<PetscReal> RecLocX() const { return mRecLocX; }
  std::vector<PetscReal> RecLocY() const { return mRecLocY; }
  std::vector<PetscReal> RecLocZ() const { return mRecLocZ; }
  std::vector<PetscReal> SrcLocX() const { return mSrcLocX; }
  std::vector<PetscReal> SrcLocY() const { return mSrcLocY; }
  std::vector<PetscReal> SrcLocZ() const { return mSrcLocZ; }
  std::vector<PetscReal> SrcRickerAmplitude() const { return mSrcRickerAmplitude; }
  std::vector<PetscReal> SrcRickerCenterFreq() const { return mSrcRickerCenterFreq; }
  std::vector<PetscReal> SrcRickerTimeDelay() const { return mSrcRickerTimeDelay; }

  std::vector<std::string> RecNames() const { return mRecNames; }
  std::vector<std::string> MovieFields() const { return mMovieFields; }

  std::vector<std::string> HomogeneousDirichlet() const { return mHomogeneousDirichletBoundaries; }

  /* Setters (mainly for testing). */
  void SetDimension(const PetscInt dim) { mNumDim = dim; }
  void SetSourceType(const std::string type) { mSourceType = type; }
  void SetReceiverFileName(const std::string type) { mReceiverFileName = type; }

};
