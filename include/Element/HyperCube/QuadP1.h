#pragma once

// stl.
#include <iostream>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>

// salvus.
#include <Utilities/Types.h>

class QuadP1 {

 private:

  const static PetscInt mNumDim = 2;
  const static PetscInt mNumVtx = 4;

 public:

  /// Empty constructor and destructor.
  QuadP1() {};
  ~QuadP1() {};

  /*
   * Linearly interpolate a value at the vertices to some location on the interior.
   * @param [in] r Reference coordinate r.
   * @param [in] s Reference coordiante s.
   * @returns A 4-D vector containing interpolation coefficients.
   */
  static RealVec4 interpolateAtPoint(const PetscReal r, const PetscReal s);

  /*
   * Check whether a given location is within, or on the edge of, an element.
   * This function will return TRUE if a point is located anywhere within an element. For
   * example, if it is located on a vertex, or edge, it will still be true. It is up to the
   * user to ensure that cases where multiple elements return true are properly handled.
   * @param [in] x Real coordinate x.
   * @param [in] y Real coordinate y.
   * @param [in] vtx Matrix containing element vertices.
   * @returns True/False based on whether or not point (x, y) is located within the element.
   */
  static bool checkHull(const PetscReal x, const PetscReal y,
                        const Eigen::Ref<const QuadVtx> &vtx);

  /**
   * Get the position of a point w.r.t. the element's reference coordiantes.
   * Uses Newton's method to do find the point.
   * @param [in] x Real coordinate x.
   * @param [in] y Real coordinate y.
   * @param [in] vtx Matrix containing element vertices.
   * @returns 2-D vector containing (x,y) mapped to (r,s).
   */
  static RealVec2 inverseCoordinateTransform(const PetscReal x, const PetscReal y,
                                             const Eigen::Ref<const QuadVtx> &vtx);

  /**
   * Get the inverse jacobian matrix, as well as the determinant of the Jacobian, at a point.
   * @param [in] r Reference coordiante r.
   * @param [in] s Reference coordinate s.
   * @param [in] vtx Matrix containing element vertices.
   * @returns A tuple (inverseJacobian, determinantJacobian).
   */
  static std::tuple<RealMat2x2, PetscReal> inverseJacobianAtPoint(
      const PetscReal r, const PetscReal s, const Eigen::Ref<const QuadVtx> &vtx,
      PetscReal &detJac, Eigen::Ref<RealMat2x2> invJac);

  /**
   * Builds the real space locations of integration points.
   * @param [in] intCrdR Integration coordiantes in r direction.
   * @param [in] intCrdS Integration coordiantes in s direction.
   * @returns A tuple of coordinates in the (x,y) realspace directions.
   */
  static std::tuple<RealVec, RealVec> buildNodalPoints(
      const Eigen::Ref<const RealVec>& intCrdR,
      const Eigen::Ref<const RealVec>& intCrdS,
      const Eigen::Ref<const QuadVtx>& vtx);

  /// Class name
  const static std::string name() { return "QUADP1"; }

};


