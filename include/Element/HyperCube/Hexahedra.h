#pragma once

// stl.
#include <map>
#include <vector>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>

// forward decl.
class Mesh;
class Source;
class Element;
class Options;
class Receiver;
class ExodusModel;

template <typename ConcreteHex>
class Hexahedra: public ConcreteHex {

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
  int mNumIntPtsT;

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
  static Eigen::MatrixXd mGrd;
  static Eigen::MatrixXd mGrdT;
  static Eigen::MatrixXd mGrdWgt;
  static Eigen::MatrixXd mGrdWgtT;
  
  // Material parameters.
  std::map<std::string,Eigen::Matrix<double,8,1>> mPar;

  // Sources and receivers.
  std::vector<std::shared_ptr<Source>> mSrc;
  std::vector<std::shared_ptr<Receiver>> mRec;

  // precomputed values for stiffness routine
  Eigen::VectorXd mDetJac;
  std::vector<Eigen::Matrix3d> mInvJac;
  
 public:

  Hexahedra<ConcreteHex>(std::unique_ptr<Options> const &options);

  /**
   * Returns the quadrature locations for a given polynomial order.
   * @param [in] order The polynmomial order.
   * @returns Vector of GLL points.
   */
  static Eigen::VectorXd GllPoints(const int order);

  /**
   * Returns the quadrature intergration weights for a polynomial order.
   * @param [in] order The polynomial order.
   * @returns Vector of quadrature weights.
   */
  static Eigen::VectorXd GllIntegrationWeights(const int order);

  /**
   * Returns the mapping from the PETSc to Salvus closure.
   * @param [in] order The polynomial order.
   * @returns Vector containing the closure mapping (field(closure(i)) = petscField(i))
   */
  static Eigen::VectorXi ClosureMapping(const int order, int elem_num,
                                        DM &distributed_mesh);

  /**
   * Setup the auto-generated gradient operator, and stores the result in mGrd.
   * @param [in] order The polynomial order.
   */
  static Eigen::MatrixXd setupGradientOperator(const int order);

  /**
   * 
   */
  static Eigen::VectorXd interpolateLagrangePolynomials(const double r,
                                                        const double s,
                                                        const double t,
                                                        const int order);
  
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
  Eigen::MatrixXd computeGradient(const Eigen::Ref<const Eigen::VectorXd>& field);

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
   * precompute constants needed for stiffness routine
   */
  void precomputeConstants();
  
  
  /**
   * Test the full stiffness routine for acoustic
   */
  Eigen::VectorXd computeStiffnessFull(const Eigen::Ref<const Eigen::VectorXd> &field, Eigen::VectorXd &mVp);
  
  /**
   * Figure out and set boundaries.
   * @param [in] mesh The mesh instance.
   */
  void setBoundaryConditions(std::unique_ptr<Mesh> const &mesh);

  /**
   * Integrate a field over the element, returning a scalar.
   * @param [in] field The field to integrate.
   */
  double integrateField(const Eigen::Ref<const Eigen::VectorXd>& field);

  /**
   * Attach the (4) vertex coordinates to the element.
   * @param [in] distributed_mesh The PETSc DM.
   */
  void attachVertexCoordinates(std::unique_ptr<Mesh> const &mesh);

  /**
   * Attach some abstract source instance to the element.
   * Test to see whether or not the source exists in the current element. If it does,
   * save a pointer to the source.
   * @param [in] sources A vector of source objects.
   */
  bool attachSource(std::unique_ptr<Source> &source, const bool finalize);

  /**
   * Attach some abstract receiver instance to the element.
   * Test to see whether or not the receiver exists in the current element. If it does,
   * save a pointer to the receiver.
   * @param [in] sources A vector of receiver objects.
   */
  void attachReceiver(std::vector<std::unique_ptr<Receiver>> receivers);

  /**
   * If an element is detected to be on a boundary, apply the Dirichlet condition to the
   * dofs on that boundary.
   * @param [in] mesh The mesh instance.
   * @param [in] options The options class.
   * @param [in] fieldname The field to which the boundary must be applied.
   */
  void applyDirichletBoundaries(std::unique_ptr<Mesh> const &mesh, std::unique_ptr<Options> const &options, const std::string &fieldname);

  /**
   *
   */
  Eigen::VectorXd getDeltaFunctionCoefficients(const double r, const double s, const double t);

  /**
   * Given a model, save the material parameters at the element vertices.
   * @param [in] model The model containing the material parameters.
   * @param [in] parameter The parameter to save.
   */
  void attachMaterialProperties(std::unique_ptr<ExodusModel> const &model, std::string parameter);

  /** Return the estimated CFL constant for the current order
   * @return The CFL estimate
   */
  double CFL_constant();
  
  /** Return the estimated element radius
   */
  virtual double estimatedElementRadius();
  
  /**
   * Given some field at the GLL points, interpolate the field to some general point.
   * @param [in] pnt Position in reference coordinates.
   */
  Eigen::MatrixXd interpolateFieldAtPoint(const Eigen::VectorXd &pnt) { return Eigen::MatrixXd(1, 1); }

  // Setters.
  inline void SetNumNew(const PetscInt num) { mElmNum = num; }
  inline void SetVtxCrd(const Eigen::Ref<const Eigen::Matrix<double,8,3>> &v) { mVtxCrd = v; }

  // Getters.
  inline bool BndElm() const { return mBndElm; }
  inline int NumDim() const { return mNumDim; }
  inline PetscInt ElmNum() const { return mElmNum; }
  inline int NumIntPnt() const { return mNumIntPnt; }
  inline int NumDofVol() const { return mNumDofVol; }
  inline int NumDofFac() const { return mNumDofFac; }
  inline int NumDofEdg() const { return mNumDofEdg; }
  inline int NumDofVtx() const { return mNumDofVtx; }
  inline Eigen::MatrixXi ClsMap() const { return mClsMap; }
  inline Eigen::MatrixXd VtxCrd() const { return mVtxCrd; }
  std::vector<std::shared_ptr<Source>> Sources() { return mSrc; }

  // Delegates.
  std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> buildNodalPoints() {
    return ConcreteHex::buildNodalPoints(mIntCrdR, mIntCrdS, mIntCrdT, mVtxCrd);
  };
  
};
