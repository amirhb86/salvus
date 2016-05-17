#pragma once

#include <petsc.h>
#include <Eigen/Dense>
template <typename ConcreteHex>
class HexahedraNew: public ConcreteHex {

/**
 * Base class of an abstract eight node hexahedra. The reference
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


private:

  const static int mNumDim = 3;
  const static int mNumVtx = 8;

  // Workspace.
  Eigen::VectorXd mDetJac;
  Eigen::VectorXd mParWork;
  Eigen::VectorXd mStiffWork;
  Eigen::MatrixXd mGradWork;

  // On Boundary.
  bool mBndElm;
  std::map<std::string,std::vector<int>> mBnd;

  // Instance variables.
  PetscInt mElmNum;
  int mPlyOrd;
  int mNumIntPnt;
  int mNumDofVtx;
  int mNumDofEdg;
  int mNumDofFac;
  int mNumDofVol;
  int mNumIntPtsR;
  int mNumIntPtsS;

  // Vertex coordinates.
  Eigen::Matrix<double,mNumVtx,mNumDim> mVtxCrd;

  // Element center.
  Eigen::Vector3d mElmCtr;

  // Closure mapping.
  Eigen::VectorXi mClsMap;

  // Quadrature parameters.
  Eigen::VectorXd mIntCrdR;
  Eigen::VectorXd mIntCrdS;
  Eigen::VectorXd mIntCrdT;
  Eigen::VectorXd mIntWgtR;
  Eigen::VectorXd mIntWgtS;
  Eigen::VectorXd mIntWgtT;

  // Matrix holding gradient information.
  Eigen::MatrixXd mGrd;

  // Material parameters.
  std::map<std::string,Eigen::Matrix<double,8,1>> mPar;

  // Sources and receivers.
  std::vector<std::shared_ptr<Source>> mSrc;
  std::vector<std::shared_ptr<Receiver>> mRec;

 public:

  HexahedraNew<ConcreteHex>(Options options);

  /**
   * Returns the quadrature locations for a given polynomial order.
   * @param [in] order The polynmomial order.
   * @returns Vector of GLL points.
   */
  static VectorXd GllPoints(const int order);

  /**
   * Returns the quadrature intergration weights for a polynomial order.
   * @param [in] order The polynomial order.
   * @returns Vector of quadrature weights.
   */
  static VectorXd GllIntegrationWeight(const int order);

  /**
   * Returns the mapping from the PETSc to Salvus closure.
   * @param [in] order The polynomial order.
   * @returns Vector containing the closure mapping (field(closure(i)) = petscField(i))
   */
  static Eigen::VectorXi ClosureMapping(const int order,
                                        DM &distributed_mesh);

  /**
   * Setup the auto-generated gradient operator, and stores the result in mGrd.
   * @param [in] order The polynomial order.
   */
  static Eigen::MatrixXd setupGradientOperator(const int order);

  /**
   * Returns an optimized stride along the r direction.
   * @param [in] f Function defined at GLL points.
   * @param [in] s_ind Index of GLL points along s-axis.
   * @param [in] numPtsS Number of integration points along the s-axis.
   * @param [in] numPtsR Number of integration points along the r-axis.
   */
  static Eigen::VectorXd rVectorStride(const Eigen::Ref<const Eigen::VectorXd>& f,
                                       const int s_ind, const int t_ind,
                                       const int numPtsR, const int numPtsS,
                                       const int numPtsT);
  /**
   * Returns an optimized stride along the s direction.
   * @param [in] f Function defined at GLL points.
   * @param [in] s_ind Index of GLL points along s-axis.
   * @param [in] numPtsS Number of integration points along the s-axis.
   * @param [in] numPtsR Number of integration points along the r-axis.
   *
   */
  static Eigen::VectorXd sVectorStride(const Eigen::Ref<const Eigen::VectorXd>& f,
                                       const int r_ind, const int t_ind,
                                       const int numPtsR, const int numPtsS,
                                       const int numPtsT);

  /**
   * Returns an optimized stride along the s direction.
   * @param [in] f Function defined at GLL points.
   * @param [in] s_ind Index of GLL points along s-axis.
   * @param [in] numPtsS Number of integration points along the s-axis.
   * @param [in] numPtsR Number of integration points along the r-axis.
   *
   */
  static Eigen::VectorXd tVectorStride(const Eigen::Ref<const Eigen::VectorXd>& f,
                                       const int r_ind, const int s_ind,
                                       const int numPtsR, const int numPtsS,
                                       const int numPtsT);

  /**
   * Compute the gradient of a field at all GLL points.
   * @param [in] field Field to take the gradient of.
   */
  Eigen::MatrixXd computeGradient(const Eigen::Ref<const Eigen::MatrixXd>& field);

  /**
   * Interpolate a parameter from vertex to GLL point.
   * @param [in] par Parameter to interpolate (i.e. VP, VS).
   */
  Eigen::VectorXd ParAtIntPts(const std::string& par);


  /**
   * Multiply a field by the test functions and integrate.
   * @param [in] f Field to calculate on.
   */
  Eigen::VectorXd applyTestAndIntegrate(const Eigen::Ref<const Eigen::VectorXd>& f);

  /**
   * Multiply a field by the gradient of the test functions and integrate.
   * @param [in] f Field to calculate on.
   */
  Eigen::VectorXd applyGradTestAndIntegrate(const Eigen::Ref<const Eigen::MatrixXd>& f);

  /**
   * Figure out and set boundaries.
   * @param [in] mesh The mesh instance.
   */
  void setBoundaryConditions(Mesh *mesh);

  /**
   * Integrate a field over the element, returning a scalar.
   * @param [in] field The field to integrate.
   */
  double integrateField(const Eigen::Ref<const Eigen::VectorXd>& field);

  /**
   * Attach the (4) vertex coordinates to the element.
   * @param [in] distributed_mesh The PETSc DM.
   */
  void attachVertexCoordinates(DM &distributed_mesh);

  /**
   * Attach some abstract source instance to the element.
   * Test to see whether or not the source exists in the current element. If it does,
   * save a pointer to the source.
   * @param [in] sources A vector of source objects.
   */
  void attachSource(std::vector<std::shared_ptr<Source>> sources);

  /**
   * Attach some abstract receiver instance to the element.
   * Test to see whether or not the receiver exists in the current element. If it does,
   * save a pointer to the receiver.
   * @param [in] sources A vector of receiver objects.
   */
  void attachReceiver(std::vector<std::shared_ptr<Receiver>> &receivers);

  /**
   * If an element is detected to be on a boundary, apply the Dirichlet condition to the
   * dofs on that boundary.
   * @param [in] mesh The mesh instance.
   * @param [in] options The options class.
   * @param [in] fieldname The field to which the boundary must be applied.
   */
  void applyDirichletBoundaries(Mesh *mesh, Options &options, const std::string &fieldname);

  /**
   * Given some field at the GLL points, interpolate the field to some general point.
   * @param [in] pnt Position in reference coordinates.
   */
  Eigen::MatrixXd interpolateFieldAtPoint(const Eigen::VectorXd &pnt) { return Eigen::MatrixXd(1, 1); }

};


