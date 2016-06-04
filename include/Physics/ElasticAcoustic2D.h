#pragma once

// stl.
#include <iostream>
#include <vector>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>

// forward decl.
class Mesh;
class Options;
class ExodusModel;

template <typename BasePhysics>
class ElasticToAcoustic2D: public BasePhysics {

 private:

  std::vector<PetscInt> mEdg, mNbr;
  std::vector<Eigen::Vector2d> mNbrCtr;

 public:

  /**** Initializers ****/
  ElasticToAcoustic2D<BasePhysics>(Options options);
  void setBoundaryConditions(Mesh *mesh);

  std::vector<std::string> PullElementalFields() const;
  Eigen::MatrixXd computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd>& u);

};
