#pragma once

// stl.
#include <iostream>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>
#include <Utilities/Types.h>

class TetP1 {

 private:

  const static int mNumDim = 3;
  const static int mNumVtx = 4;
  
 public:

  TetP1() {};
  static Eigen::Vector4d interpolateAtPoint(double r, double s, double t);

  /**
   * Checks whether a given point in realspace (x, z) is within the current element.
   * A simple convex hull algorithm is implemented. Should work as long as the sides of the element are straight
   * lines, but will likely fail for higher order shape functions.
   * @param [in] x X-coordinate in real space.
   * @param [in] z Z-coordinate in real space.
   * @returns True if point is inside, False if not.
   */
  static bool checkHull(double x, double y, double z,
                        const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &vtx);
  
  /**
   * Given a point in realspace, determines the equivalent location in the reference element.
   * Since the shape function is linear, the inverse transform is a simple analytic formula.
   * @param [in] x_real X-coordinate in real space.
   * @param [in] y_real Y-coordinate in real space.
   * @param [in] z_real Z-coordinate in real space.
   * @param [in] vtx Coordinates of the 4 vertices.
   * @return A Vector (r,s,t) containing the coordinates in the reference element.
   */
  static Eigen::Vector3d inverseCoordinateTransform(double x, double y, double z,
                                                    const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &vtx);

  /**
   * 3x3 inverse Jacobian matrix.  This method returns an Eigen::Matrix
   * representation of the inverse Jacobian for a given element with vtx
   * @param [in] vtx Vertices of element
   * @returns (inverse Jacobian matrix,determinant of that matrix) as a `std::tuple`. Tuples can be
   * "destructured" using a `std::tie`.
   */
  static std::tuple<Eigen::Matrix3d, PetscReal> inverseJacobian(
                                                                const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &vtx);

  static std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> buildNodalPoints(
      const Eigen::Ref<const Eigen::VectorXd>& intCrdR,
      const Eigen::Ref<const Eigen::VectorXd>& intCrdS,
      const Eigen::Ref<const Eigen::VectorXd>& intCrdT,
      const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>>& vtx);


  /// Class name
  const static ElementType type() { return ElementType::TETP1; }
  
};


