#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <Utilities/Types.h>

class Mesh;
class Options;

template <typename Base>
class HomogeneousDirichlet: public Base {

 private:

  std::vector<PetscInt> mBndVtx;
  std::vector<PetscInt> mBndEdg;
  std::vector<PetscInt> mBndFac;

 public:

  HomogeneousDirichlet<Base>(std::unique_ptr<Options> const &options);

  void setBoundaryConditions(std::unique_ptr<Mesh> const &mesh);

  RealMat computeStiffnessTerm(const Eigen::Ref<const RealMat>& u);

  const static std::string Name() { return "HomogeneousDirichlet_" + Base::Name(); }

};


