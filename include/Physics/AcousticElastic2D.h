#pragma once

// std.
#include <iostream>
#include <vector>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>

class Mesh;
class Options;

template <typename BasePhysics>
class AcousticToElastic2D: public BasePhysics {

 private:

  const static int num_scalar_components = 1;
  const static int num_elastic_components = 2;
  std::vector<PetscInt> mEdg;
  Eigen::Matrix<double,Eigen::Dynamic,num_elastic_components> uElastic;
  Eigen::Matrix<double,Eigen::Dynamic,num_scalar_components> uScalar;

 public:

  /**** Initializers ****/
  AcousticToElastic2D<BasePhysics>(Options options);
  void setBoundaryConditions(Mesh *mesh);

  std::vector<std::string> PullElementalFields() const;
  Eigen::MatrixXd computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd>& u);

};