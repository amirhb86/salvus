#pragma once

// stl.
#include <iostream>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>
#include <Utilities/Types.h>

class TriP1 {

 private:

  const static int mNumDim = 2;
  const static int mNumVtx = 3;
  
 public:

  /// Empty constructor and destructor.
  TriP1() {};
  ~TriP1() {};
  
  static Eigen::Vector3d interpolateAtPoint(const double r, const double s);

  /**
   * Checks whether a given point in realspace (x, z) is within the current element.
   * A simple convex hull algorithm is implemented. Should work as long as the sides of the element are straight
   * lines, but will likely fail for higher order shape functions.
   * @param [in] x X-coordinate in real space.
   * @param [in] z Z-coordinate in real space.
   * @returns True if point is inside, False if not.
   */
  static bool checkHull(const double x, const double y,
                        const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &vtx);
  
  /**
   * Given a point in realspace, determines the equivalent location in the reference element.
   * Since the shape function is linear, the inverse transform is a simple analytic formula.
   * @param [in] x_real X-coordinate in real space.
   * @param [in] z_real Z-coordinate in real space.
   * @param [in] vtx Coordinates of the 3 vertices.
   * @return A Vector (eps, eta) containing the coordinates in the reference element.
   */
  static Eigen::Vector2d inverseCoordinateTransform(const double x, const double y,
                                                    const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &vtx);

  /**
   * 2x2 inverse Jacobian matrix.  This method returns an Eigen::Matrix
   * representation of the inverse Jacobian for a given element with vtx
   * @param [in] vtx Vertices of element
   * @returns (inverse Jacobian matrix,determinant of that matrix) as a `std::tuple`. Tuples can be
   * "destructured" using a `std::tie`.
   */
  static std::tuple<Eigen::Matrix2d, PetscReal> inverseJacobian(
      const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &vtx);

  static std::tuple<Eigen::VectorXd, Eigen::VectorXd> buildNodalPoints(
      const Eigen::Ref<const Eigen::VectorXd>& intCrdR,
      const Eigen::Ref<const Eigen::VectorXd>& intCrdS,
      const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>>& vtx);

  /// Class name
  const static std::string Name() { return "TriP1"; }

  /// Class type
  const static ElementType type() { return ElementType::TRIP1; }
  
};


