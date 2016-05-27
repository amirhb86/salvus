#pragma once

// stl.
#include <iostream>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>

class QuadP1 {

 private:

  const static int mNumDim = 2;
  const static int mNumVtx = 4;

 public:

  QuadP1() {};
  static Eigen::Vector4d interpolateAtPoint(const double r, const double s);
  static bool checkHull(const double x, const double y,
                        const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &vtx);
  static Eigen::Vector2d inverseCoordinateTransform(const double x, const double y,
                                                    const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &vtx);
  static std::tuple<Eigen::Matrix2d, PetscReal> inverseJacobianAtPoint(
      const PetscReal r, const PetscReal s,
      const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &vtx);
  static std::tuple<Eigen::VectorXd, Eigen::VectorXd> buildNodalPoints(
      const Eigen::Ref<const Eigen::VectorXd>& intCrdR,
      const Eigen::Ref<const Eigen::VectorXd>& intCrdS,
      const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>>& vtx);
};


