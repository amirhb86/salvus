#pragma once

#include <Eigen/Dense>
#include <Model/ExodusModel.h>
#include <Source/Source.h>
#include <Mesh/Mesh.h>
#include <Element/Element.h>

extern "C" {
#include "Quad/Autogen/quad_autogen.h"
};

using namespace Eigen;

/**
 * Base class of an abstract four node quadrilateral. The reference
 * element is set up as below, where vertices are ordered with
 * right-hand-rule normals pointing out.
 *
 *                                  (v7)                    (v6)
 *                                    /---------------------/+
 *                                 /--|                  /-- |
 *                             /---   |              /---    |
 *                          /--       |           /--        |
 *                       /--          |         --           |
 *                (v4) +------------------------+ (v5)       |
 *                     |              |         |            |
 *                     |              |         |            |
 *                     |              |         |            |
 *   ^                 |              |         |            |
 *   | (t)             |              |         |            |
 *   |                 |            --+-------- |------------+ (v2)
 *   |                 |        ---/ (v1)       |         /--
 *   |                 |     --/                |     /---
 *   |         (s)     | ---/                   |  /--
 *   |         /-      +/-----------------------+--
 *   |     /---       (v0)                    (v3)
 *   | /---
 *   +-----------------> (r)
 *
 */
class Hexahedra: public Element {

  /****************************************************************************
   *                         P1 Shape Functions.
   *****************************************************************************/

  /**
   * Checks whether a given point in realspace (x, z) is within the current element.
   * A simple convex hull algorithm is implemented. Should work as long as the sides of the element are straight
   * lines, but will likely fail for higher order shape functions.
   * @param [in] x X-coordinate in real space.
   * @param [in] z Z-coordinate in real space.
   * @returns True if point is inside, False if not.
   */
  bool mCheckHull(double x, double z);  

protected:

  int mNumIntPtsR;
  /** < Num of integration points in the `r` direction (e.g. 5 for a 4th order gll basis) */
  int mNumIntPtsS;
  /** < Num of integration points in the `s` direction (e.g. 5 for a 4th order gll basis) */
  int mNumIntPtsT;
  /** < Num of integration points in the `t` direction (e.g. 5 for a 4th order gll basis) */
  MatrixXd mGrd;
  /** < Derivative of shape function n (col) at pos. m (row) */
  VectorXd mIntWeightsR;
  /** < Integration weights along `r` direction. */
  VectorXd mIntWeightsS;
  /** < Integration weights along `s` direction. */
  VectorXd mIntWeightsT;
  /** < Integration weights along `t` direction. */
  VectorXd mIntCoordR; /** < Integration points along `r` direction */
  VectorXd mIntCoordS;  /** < Integration points along `s` direction */
  VectorXd mIntCoordT;/** < Integration points along `t` direction */

public:

  /**
   * TODO: Write...
   *
   * @param [in] function Field defined across an element, for which a view is desired.
   * @param [in] r_index Eta_index in the reference element
   * @returns A const pointer with the proper stride and starting point for the desired field points.
   */
  Eigen::Map<const Eigen::VectorXd> rVectorStride(const Eigen::VectorXd &function,
                                                  const int &s_index,
                                                  const int &t_index);

  /**
   * TODO: Write...
   *
   * @param [in] function Field defined across an element, for which a view is desired.
   * @param [in] s_index Eps_index in the reference element
   * @returns A const pointer with the proper stride and starting point for the desired field points.
   */
  Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>> sVectorStride(const Eigen::VectorXd &function,
                                                                           const int &r_index,
                                                                           const int &t_index);

  /**
   * * TODO: Write...
   * @param [in] function Field defined across an element, for which a view is desired.
   * @param [in] t_index Eps_index in the reference element
   * @returns A const pointer with the proper stride and starting point for the desired field points.
   */
  Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>> tVectorStride(const Eigen::VectorXd &function,
                                                                           const int &r_index,
                                                                           const int &s_index);

  // TODO: HEXP1
  /**
   * 3x3 inverse Jacobian matrix at a point (r, s, t).  This method returns an Eigen::Matrix
   * representation of the inverse Jacobian at a particular point.
   * @param [in] r `r` position on the reference element.
   * @param [in] s `s` position on the reference element.
   * @param [in] t `t` position on the reference element.
   * @returns (inverse Jacobian matrix,determinant of that matrix) as a `std::tuple`. Tuples can be
   * "destructured" via `std::tie(t1,t2) = inverseJacobian()`.
   */
  std::tuple<Matrix3d, PetscReal> inverseJacobianAtPoint(PetscReal r, PetscReal s, PetscReal t);

  // TODO: HEXP1
  /**
   * Given a point in realspace, determines the equivalent location in the reference element.
   * Since the shape function are bilinear, the cross terms introduce nonlinearities into the shape function
   * interpolation. As such, we used a simple implementation of Newton's method to minimize the objective function:
   * q = x_real - shape_functions * vertex_coordinates(r,s,t), for each coordinate (x,y,z).
   * @param [in] x_real X-coordinate in real space.
   * @param [in] y_real Y-coordinate in real space.
   * @param [in] z_real Z-coordinate in real space.
   * @return A Vector (r,s,t) containing the coordinates in the reference element.
   */
  Vector3d inverseCoordinateTransform(const double &x_real, const double &y_real, const double &z_real);
  
  /**
   * Constructor.
   * Sets quantities such as number of dofs, among other things, from the options class.
   * @param [in] options Populated options class.
   */
  Hexahedra(Options options);

  /****************************************************************************
   *                       STATIC UTILITY FUNCTIONS
   *****************************************************************************/

  /**
   * Returns the quadrature locations for a given polynomial order.
   * @param [in] order The polynmomial order.
   * @returns Vector of GLL points.
   */
  static VectorXd GllPointsForOrder(const int order);

  /**
   * Returns the quadrature intergration weights for a polynomial order.
   * @param [in] order The polynomial order.
   * @returns Vector of quadrature weights.
   */
  static VectorXd GllIntegrationWeightForOrder(const int order);
  
  /**
   * Provides interpolation vector for quadrilaterals with vertex velocities based on right hand rule.
   * @param [in] r `r` in the reference element.
   * @param [in] s `s` in the reference element.
   * @param [in] t `t` in the reference element.
   * @returns A vector containing the [4] coefficients from each shape function.
   * Usage: double velocity_at_r_s_t = interpolateAtPoint(r,s,t).dot(mMaterialVelocityAtVertices);
   */
  static VectorXd interpolateAtPoint(double r, double s, double t);

  /**
   * Attaches a material parameter to the vertices on the current element.
   * Given an exodus model object, we use a kD-tree to find the closest parameter to a vertex. In practice, this
   * closest parameter is often exactly coincident with the vertex, as we use the same exodus model for our mesh
   * as we do for our parameters.
   * @param [in] model An exodus model object.
   * @param [in] parameter_name The name of the field to be added (i.e. velocity, c11).
   * @returns A Vector with 4-entries... one for each Element vertex, in the ordering described above.
   */
  VectorXd __attachMaterialProperties(ExodusModel *model,
                                      std::string parameter_name);

  /**
   * Utility function to integrate a field over the element. This could probably be made static, but for now I'm
   * just using it to check some values.
   * @param [in] field The field which to integrate, defined on each of the gll points.
   * @returns The scalar value of the field integrated over the element.
   */
  double integrateField(const Eigen::VectorXd &field);

  /**
   * Setup the auto-generated gradient operator, and stores the result in mGrd.
   */
  void setupGradientOperator();


  /**
   * Queries the passed DM for the vertex coordinates of the specific element. These coordinates are saved
   * in mVertexCoordiantes.
   * @param [in] distributed_mesh PETSc DM object.
   */
  void attachVertexCoordinates(DM &distributed_mesh);

  /**
   * Attach source.
   * Given a vector of abstract source objects, this function will query each for its spatial location. After
   * performing a convex hull test, it will perform a quick inverse problem to determine the position of any sources
   * within each element in reference coordinates. These reference coordinates are then saved in the source object.
   * References to any sources which lie within the element are saved in the mSrc vector.
   * @param [in] sources A vector of all the sources defined for a simulation run.
   */
  void attachSource(std::vector<std::shared_ptr<Source>> sources);

  /**
   * Returns the mapping from the PETSc to Salvus closure.
   * @param [in] order The polynomial order.
   * @param [in] dimension Element dimension.
   * @returns Vector containing the closure mapping (field(closure(i)) = petscField(i))
   */
  VectorXi ClosureMapping(const int order, const int dimension, DM &distributed_mesh);

  /**
   * Builds element's closure mapping
   * @param [in] mesh The mesh object needed to get surface and edge orientations
   */
  void BuildClosureMapping(DM &distributed_mesh);
  
  /**
   * Builds nodal coordinates (x,y,z) on all mesh degrees of freedom.
   */
  std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> buildNodalPoints();

  void applyDirichletBoundaries(Mesh *mesh, Options &options, const std::string &fieldname);
  
  // for testing
  inline int __GetNumIntPtsR() { return mNumIntPtsR; }
  inline int __GetNumIntPtsS() { return mNumIntPtsS; }
  inline int __GetNumIntPtsT() { return mNumIntPtsT; }
  inline VectorXd __GetIntCoordR() { return mIntCoordR; }
  inline VectorXd __GetIntCoordS() { return mIntCoordS; }
  inline VectorXd __GetIntCoordT() { return mIntCoordT; }
  
  // to get it to compile...
  void prepareStiffness();
  void assembleElementMassMatrix(Mesh *mesh);
  MatrixXd computeSourceTerm(double time);
  MatrixXd computeStiffnessTerm(const MatrixXd &displacement);
  std::vector<std::string> PullElementalFields() const { return {"u"}; }
  std::vector<std::string> PushElementalFields() const { return {"a"}; }
  std::vector<std::string> PushElementalFields();  
  void attachMaterialProperties(ExodusModel *model);
  
};
