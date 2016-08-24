#pragma once
#include <hdf5.h>
#include <iostream>
#include <Eigen/Dense>
#include <Utilities/Types.h>

template <typename Element>
class SaveSurfaceGreenFunction: public Element {

  static MPI_Group mMpiGroup;
  static MPI_Comm mMpiComm;

 public:

  SaveSurfaceGreenFunction<Element>(std::unique_ptr<Options> const &options);
  void recordDynamicFields(const Eigen::Ref<const RealMat>& field);

  const static std::string Name() { return "SaveSurface_" + Element::Name(); };

};

