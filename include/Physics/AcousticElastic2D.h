#pragma once

// std.
#include <iostream>
#include <vector>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>

class Mesh;
class Options;
class ExodusModel;

template <typename BasePhysics>
class AcousticToElastic2D: public BasePhysics {

 private:

  std::vector<double> mRho_0;
  std::vector<PetscInt> mEdg, mNbr;
  std::vector<Eigen::Vector2d> mNbrCtr;

 public:

  /**** Initializers ****/
  AcousticToElastic2D<BasePhysics>(std::unique_ptr<Options> const &options);
  void setBoundaryConditions(Mesh *mesh);

  void attachMaterialProperties(std::unique_ptr<ExodusModel> const &model);
  std::vector<std::string> PullElementalFields() const;
  Eigen::MatrixXd computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd>& u);

};