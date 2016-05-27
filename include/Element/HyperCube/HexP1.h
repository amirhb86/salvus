#pragma once

// stl.
#include <iostream>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>

class HexP1 {

 private:

  const static int mNumDim = 3;
  const static int mNumVtx = 8;

 public:

  HexP1() {};
  static Eigen::VectorXd interpolateAtPoint(const double r, const double s, const double t);
  static bool checkHull(const double x, const double y, const double z,
                        const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &vtx);
  static Eigen::Vector3d inverseCoordinateTransform(const double x, const double y, const double z,
                                                    const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &vtx);
  static std::tuple<Eigen::Matrix3d, PetscReal>
  inverseJacobianAtPoint(const PetscReal r, const PetscReal s, const PetscReal t,                         
                         const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &vtx);
  
  static std::tuple<Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd>
  buildNodalPoints(const Eigen::Ref<const Eigen::VectorXd>& intCrdR,
                   const Eigen::Ref<const Eigen::VectorXd>& intCrdS,
                   const Eigen::Ref<const Eigen::VectorXd>& intCrdT,
                   const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>>& vtx);
  
};


