#pragma once

// stl.
#include <map>
#include <memory>
#include <vector>
#include <iostream>
#include <stdexcept>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>

// salvus
#include <Utilities/Types.h>

// forward decl.
class Mesh;
class Source;
class Options;
class Receiver;
class ExodusModel;

template <typename ConcreteShape>
class TensorQuad: public ConcreteShape {

  /** \class TensorQuad
   *
   * \brief Base class for all Quads with a tensorized basis.
   *
   * This class is at the base mixin level for all tensorized TensorQuad elements. It should concretely define things common
   * to all Quads. For example, here we define specific methods to attach vertex coordinates to each element, as well
   * as sources and receivers. These last two methods bring up an interesting point however. As you might notice, the
   * mixin approach lets us be perform class nesting to an arbitrary depth, without incurring any runtime cost. So,
   * for example, we could have a super class called something like Element2D, where source and receiver attachments
   * are implemented. However, at one point we need to decide on an optimal level of complexity. Since TensorQuad is the
   * absolute base for all further derivations, this means that attachReceiver() and attachSource() will need to be
   * duplicated for Triangles, and for all other 2D elements. Is this an acceptable level of code duplication, in order
   * to cut down on the class heiarchy depth? Perhaps. Either way, we should discuss this, and also note that following
   * these mixin principles, things are easily modified.
   *
   * To initialize a new TensorQuad element, you'll need at least one template parameter. For instance:
   *
   * \code
   * auto elem = TensorQuad<QuadP1>
   * auto elem = TensorQuad<Acoustic<QuadP1>>
   * auto elem = TensorQuad<Gravity<Attenuation<Acoustic<QuadP1>>>>
   * \endcode
   *
   * are all valid (provided that the relevant physics classes are written). Just as an example of what it would take
   * for _full_ generality, as discussed above, we might have:
   *
   * \code
   * auto elem = Elem2D<TensorQuad<GllQuadrature<Gravity<Attenuation<Acoustic<QuadP1>>>>>>.
   * \endcode
   *
   * Pretty insane, right? But, it would allow us to do some neat things (like fixed-size arrays for all quadrature
   * orders). So, I dunno man.
   *
   * Unless otherwise mentioned, all instances of `field` in function signatures refer to some function defiend
   * at the GLL points.
   */

private:

  // Static variables.
  const static PetscInt mNumDim = 2;
  const static PetscInt mNumVtx = 4;
  const static PetscInt mMaxOrder = 10;

  // Workspace.
  RealVec mDetJac;
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

  // Vertex coordinates.
  QuadVtx mVtxCrd;

  // Element center.
  RealVec2 mElmCtr;

  // Closure mapping.
  IntVec mClsMap;
  std::vector<PetscInt> mEdgMap;

  // Quadrature parameters.
  RealVec mIntCrdR;
  RealVec mIntCrdS;
  RealVec mIntWgtR;
  RealVec mIntWgtS;

  // Matrix holding gradient information.
  RealMat mGrd;

  // Material parameters.
  std::map<std::string,RealVec4> mPar;

  // Sources and receivers.
  std::vector<std::unique_ptr<Source>> mSrc;
  std::vector<std::unique_ptr<Receiver>> mRec;

 public:

  /// Allocates memory for work arrays, most private variables.
  TensorQuad<ConcreteShape>(std::unique_ptr<Options> const &options);

  /// All memory should be properly managed already.
  ~TensorQuad<ConcreteShape>() {};

  /**
   * Returns GLL point locations for a given polynomial order.
   * @param [in] order Polynomial order.
   * @returns Vector of GLL point locations.
   */
  static RealVec GllPointsForOrder(const PetscInt order);

  /**
   * Returns GLL integration weights for a given polynomial order.
   * @param [in] order Polynomial order.
   * @returns Vector of GLL integration weights.
   */
  static RealVec GllIntegrationWeightsForOrder(const PetscInt order);

  /**
   * Returns a matrix of the test function derivatives for a given polynomial order.
   * @param [in] order Polynomial order.
   * @returns Matrix of derivatives with the ordering described above.
   */
  static RealMat setupGradientOperator(const PetscInt order);

  /**
   * Given reference coordinates, interpolate Lagrange polynomials at a point.
   * @param [in] r Reference coordinate r.
   * @param [in] s Reference coordinate s.
   * @param [in] order Polynomial order.
   */
  static RealVec interpolateLagrangePolynomials(const PetscReal r, const PetscReal s,
                                                const PetscInt order);

  static IntVec  ClosureMappingForOrder(const PetscInt order);

  /**
   * Compute the gradient of a field at all GLL points.
   * @param [in] field Field to take the gradient of.
   * @returns nGll x nDim matrix containing field gradient components.
   */
  RealMat computeGradient(const Eigen::Ref<const RealVec>& field);

  /**
   * Interpolate a parameter from vertex to GLL point.
   * @param [in] par Parameter to interpolate (i.e. VP, VS).
   * @returns ngll vector containing material parameters.
   */
  RealVec ParAtIntPts(const std::string& par);

  /**
   * Multiply a field by the test functions and integrate.
   * @param [in] f Field to calculate on.
   * @returns Coefficients at gll points.
   */
  RealVec applyTestAndIntegrate(const Eigen::Ref<const RealVec>& f);

  /**
   * Multiply a field by the gradient of the test functions and integrate.
   * @param [in] f Field to calculate on.
   * @returns Coefficients at gll points.
   */
  RealVec applyGradTestAndIntegrate(const Eigen::Ref<const RealMat>& f);


  /**
   * Multiply a field by the test functions on a certain edge, and integrate.
   * @param [in] field Working field, defined at all GLL points.
   * @param [in] edg Edge number to integrate.
   * @returns Coefficients at all gll points (i.e. zeroes in the interior).
   */
  RealVec applyTestAndIntegrateEdge(const Eigen::Ref<const RealVec>& f,
                                    const PetscInt edg);

  /**
   * Given an edge, return a perpendicular normal vector to that edge.
   * @param [in] edg Edge number.
   * @returns A 2d "normal" pointing outwards.
   */
  Eigen::Vector2d getEdgeNormal(const PetscInt edg);

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
  bool attachReceiver(std::unique_ptr<Receiver> &receiver, const bool finalize);

  /**
   * Given a delta function at some location (r,s), computes the coefficients at the
   * GLL points that would integrate to the delta function over the element.
   * @param [in] r Reference coordinate.
   * @param [in] s Reference coordinate.
   */
  Eigen::VectorXd getDeltaFunctionCoefficients(const double r, const double s);

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
  Eigen::MatrixXd interpolateFieldAtPoint(const Eigen::VectorXd &pnt) { return Eigen::MatrixXd(1, 1); }

  // Setters.
  inline void SetNumNew(const PetscInt num) { mElmNum = num; }
  inline void SetVtxCrd(const Eigen::Ref<const Eigen::Matrix<double,4,2>> &v) { mVtxCrd = v; }
  inline void SetCplEdg(const std::vector<PetscInt> &v) { mEdgMap = v; }

  // Getters.
  inline bool BndElm()        const { return mBndElm; }
  inline PetscInt NumDim()    const { return mNumDim; }
  inline PetscInt ElmNum()    const { return mElmNum; }
  inline PetscInt NumVtx()    const { return mNumVtx; }
  inline PetscInt NumIntPnt() const { return mNumIntPnt; }
  inline PetscInt NumDofVol() const { return mNumDofVol; }
  inline PetscInt NumDofFac() const { return mNumDofFac; }
  inline PetscInt NumDofEdg() const { return mNumDofEdg; }
  inline PetscInt NumDofVtx() const { return mNumDofVtx; }
  inline IntVec ClsMap()      const { return mClsMap; }
  inline QuadVtx VtxCrd()     const { return mVtxCrd; }
  const inline std::vector<std::unique_ptr<Source>> &Sources() const { return mSrc; }

  inline static PetscInt MaxOrder() { return mMaxOrder; }

  // Delegates.
  std::tuple<Eigen::VectorXd, Eigen::VectorXd> buildNodalPoints() {
    return ConcreteShape::buildNodalPoints(mIntCrdR, mIntCrdS, mVtxCrd);
  };


  // TODO: DO WE STILL NEED THESE?
  /**
   * Figure out and set boundaries.
   * @param [in] mesh The mesh instance.
   */
  void setBoundaryConditions(std::unique_ptr<Mesh> const &mesh) {};

  /**
   * Integrate a field over the element, returning a scalar.
   * @param [in] field The field to integrate.
   */
  double integrateField(const Eigen::Ref<const Eigen::VectorXd>& field) {};

  virtual double CFL_constant() { return 1; };

  /** Return the estimated element radius
   * @return The CFL estimate
   */
  // TODO: MODIFY
  double estimatedElementRadius() { return 1; };

  /**
   * If an element is detected to be on a boundary, apply the Dirichlet condition to the
   * dofs on that boundary.
   * @param [in] mesh The mesh instance.
   * @param [in] options The options class.
   * @param [in] fieldname The field to which the boundary must be applied.
   */
  void applyDirichletBoundaries(std::unique_ptr<Mesh> const &mesh,
                                std::unique_ptr<Options> const &options,
                                const std::string &fieldname) {};



};


