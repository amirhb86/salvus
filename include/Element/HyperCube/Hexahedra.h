#pragma once

// stl.
#include <map>
#include <vector>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>
#include <Utilities/Types.h>

// Maximum order for Hexahedra. Order 8 and 9 have generated code, but
// take quite a long time to compile, which is quite annoying for
// day-to-day use and development. Compile with -DHEX_MAX_ORDER 9 to
// set otherwise.
#ifndef HEX_MAX_ORDER
#define HEX_MAX_ORDER 7
#endif


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

  const static PetscInt mNumDim = 3;
  const static PetscInt mNumVtx = 8;

  const static PetscInt mMaxOrder = HEX_MAX_ORDER; // defined above

  // Workspace.  
  RealVec mParWork;
  RealVec mStiffWork;
  RealMat mGradWork;

  // On Boundary.
  bool mBndElm;
  std::map<std::string,std::vector<PetscInt>> mBnd;

  // Instance variables.
  PetscInt mElmNum;
  PetscInt mPlyOrd;
  PetscInt mNumIntPnt;
  PetscInt mNumDofVtx;
  PetscInt mNumDofEdg;
  PetscInt mNumDofFac;
  PetscInt mNumDofVol;
  PetscInt mNumIntPtsR;
  PetscInt mNumIntPtsS;
  PetscInt mNumIntPtsT;

  // Vertex coordinates.
  HexVtx mVtxCrd;

  // Element center.
  Eigen::Vector3d mElmCtr;

  // Closure mapping.
  IntVec mClsMap;

  // Quadrature parameters.
  RealVec mIntCrdR;
  RealVec mIntCrdS;
  RealVec mIntCrdT;
  RealVec mIntWgtR;
  RealVec mIntWgtS;
  RealVec mIntWgtT;

  // Matrix holding gradient information.
  RealMat mGrd;
  RealMat mGrdT;
  RealMat mGrdWgt;
  RealMat mGrdWgtT;
  
  // Material parameters.
  std::map<std::string,RealVec> mPar;

  // Sources and receivers.
  std::vector<std::unique_ptr<Source>> mSrc;
  std::vector<std::unique_ptr<Receiver>> mRec;

  // precomputed values for stiffness routine
  RealVec mDetJac;
  std::vector<RealMat3x3> mInvJac;
  
 public:

  Hexahedra<ConcreteHex>(std::unique_ptr<Options> const &options);

  /**
   * Returns the quadrature locations for a given polynomial order.
   * @param [in] order The polynmomial order.
   * @returns Vector of GLL points.
   */
  static RealVec GllPointsForOrder(const PetscInt order);

  /**
   * Returns the quadrature intergration weights for a polynomial order.
   * @param [in] order The polynomial order.
   * @returns Vector of quadrature weights.
   */
  static RealVec GllIntegrationWeights(const PetscInt order);

  /**
   * Returns the mapping from the PETSc to Salvus closure.
   * @param [in] order The polynomial order.
   * @returns Vector containing the closure mapping (field(closure(i)) = petscField(i))
   */
  static IntVec ClosureMapping(const PetscInt order, PetscInt elem_num,
                                        DM &distributed_mesh);

  /**
   * Setup the auto-generated gradient operator, and stores the result in mGrd.
   * @param [in] order The polynomial order.
   */
  static RealMat setupGradientOperator(const PetscInt order);

  /**
   * 
   */
  static RealVec interpolateLagrangePolynomials(const PetscReal r,
                                                        const PetscReal s,
                                                        const PetscReal t,
                                                        const PetscInt order);

  /**
   * Compute the gradient of a field at all GLL points.
   * @param [in] field Field to take the gradient of.
   */
  RealMat computeGradient(const Eigen::Ref<const RealVec>& field);

  /**
   * Interpolate a parameter from vertex to GLL point.
   * @param [in] par Parameter to interpolate (i.e. VP, VS).
   */
  RealVec ParAtIntPts(const std::string& par);


  /**
   * Multiply a field by the test functions and integrate.
   * @param [in] f Field to calculate on.
   */
  RealVec applyTestAndIntegrate(const Eigen::Ref<const RealVec>& f);

  /**
   * Multiply a field by the test functions on a certain edge, and integrate.
   * @param [in] field Working field, defined at all GLL points.
   * @param [in] edg Edge number to integrate.
   * @returns Coefficients at all gll points (i.e. zeroes in the interior).
   */
  RealVec applyTestAndIntegrateEdge(const Eigen::Ref<const RealVec> &f,
                                    const PetscInt edg);

  /**
   * Multiply a field by the gradient of the test functions and integrate.
   * @param [in] f Field to calculate on.
   */
  RealVec applyGradTestAndIntegrate(const Eigen::Ref<const Eigen::MatrixXd>& f);

  /**
   * precompute constants needed for stiffness routine
   */
  void precomputeConstants();

  void setEdgeToValue(const PetscInt edg, const PetscReal val, Eigen::Ref<RealVec> f);
  
  
  /**
   * Test the full stiffness routine for acoustic
   */
  RealVec computeStiffnessFull(const Eigen::Ref<const RealVec> &field, RealVec &mVp);
  
  /**
   * Figure out and set boundaries.
   * @param [in] mesh The mesh instance.
   */
  void setBoundaryConditions(std::unique_ptr<Mesh> const &mesh) {};

  /**
   * Integrate a field over the element, returning a scalar.
   * @param [in] field The field to integrate.
   */
  PetscReal integrateField(const Eigen::Ref<const RealVec>& field);

  /**
   * Attach the (4) vertex coordinates to the element.
   * @param [in] distributed_mesh The PETSc DM.
   */
  void attachVertexCoordinates(std::unique_ptr<Mesh> const &mesh);

  /** Precompute any terms needed on the element level, e.g.,
      jacobians or velocities at nodes.
  */
  void precomputeElementTerms() {}
  
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
  bool attachReceiver(std::unique_ptr<Receiver> &receiver, const bool finalize);

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
  RealVec getDeltaFunctionCoefficients(const Eigen::Ref<RealVec>& pnt);

  /**
   * Given a model, save the material parameters at the element vertices.
   * @param [in] model The model containing the material parameters.
   * @param [in] parameter The parameter to save.
   */
  void attachMaterialProperties(std::unique_ptr<ExodusModel> const &model, std::string parameter);

  /**
   * Given some field at the GLL points, interpolate the field to some general point.
   * @param [in] pnt Position in reference coordinates.
   */
  RealMat interpolateFieldAtPoint(const RealVec &pnt) { return Eigen::MatrixXd(1, 1); }

  // Setters.
  inline void SetNumNew(const PetscInt num) { mElmNum = num; }
  inline void SetVtxCrd(const Eigen::Ref<const HexVtx> &v) { mVtxCrd = v; }

  // Getters.
  inline bool BndElm() const { return mBndElm; }
  inline PetscInt NumDim() const { return mNumDim; }
  inline PetscInt ElmNum() const { return mElmNum; }
  inline PetscInt NumIntPnt() const { return mNumIntPnt; }
  inline PetscInt NumDofVol() const { return mNumDofVol; }
  inline PetscInt NumDofFac() const { return mNumDofFac; }
  inline PetscInt NumDofEdg() const { return mNumDofEdg; }
  inline PetscInt NumDofVtx() const { return mNumDofVtx; }
  inline IntVec ClsMap() const { return mClsMap; }
  inline int PlyOrd()      const { return mPlyOrd; }
  inline RealMat VtxCrd() const { return mVtxCrd; }
  inline static PetscInt MaxOrder() { return mMaxOrder; }
  inline void SetVtxPar(const Eigen::Ref<const RealVec> &v, const std::string &par) { mPar[par] = v; }
  const inline std::vector<std::unique_ptr<Source>> &Sources() const { return mSrc; }
  const inline std::vector<std::unique_ptr<Receiver>> &Receivers() const { return mRec; }

  // Delegates.
  std::tuple<RealVec, RealVec, RealVec> buildNodalPoints() {
    return ConcreteHex::buildNodalPoints(mIntCrdR, mIntCrdS, mIntCrdT, mVtxCrd);
  };

  const static std::string Name() { return "TensorHex_" + ConcreteHex::Name(); }

};
