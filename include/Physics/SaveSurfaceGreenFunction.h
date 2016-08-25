#pragma once
#include <hdf5.h>
#include <hdf5_hl.h>
#include <iostream>
#include <Eigen/Dense>
#include <Utilities/Types.h>

template <typename Element>
class SaveSurfaceGreenFunction: public Element {

  static hid_t mFileId;
  static MPI_Comm mMpiComm;
  static MPI_Group mMpiGroup;
  static PetscInt mNumMixins;
  static PetscInt mNumTimeStep;
  static UncompressedWavefieldContainer<PetscReal> container;

  PetscInt mMixNum;
  PetscInt mNumComponents;
  PetscInt mCountTimeSteps;

  void write();

 public:

  SaveSurfaceGreenFunction<Element>(std::unique_ptr<Options> const &options);
  ~SaveSurfaceGreenFunction<Element>();
  void recordDynamicFields(const Eigen::Ref<const RealMat>& field);

  const static std::string Name() { return "SaveSurface_" + Element::Name(); };

};

class SaveSurfaceGreenFunctionTestStub {

 public:
  SaveSurfaceGreenFunctionTestStub(std::unique_ptr<Options> const &options) {};
  ~SaveSurfaceGreenFunctionTestStub() {};
  const static std::string Name() { return "TestStub"; }
  const static std::vector<std::string> PullElementalFields() { return {"ux_test", "uy_test"}; };
  static PetscInt NumIntPnt() { return 25; }

};
