#pragma once

#include <map>
#include <memory>
#include <vector>

#include <petsc.h>
#include <Eigen/Dense>

// forward decl.
class Mesh;
class Source;
class Options;
class Receiver;
class ExodusModel;

template <typename ConcreteShape>
class Quad: public ConcreteShape {

  /** \class Quad
   *
   * \brief Base class for all Quads with a tensorized basis.
   *
   * This class is at the base mixin level for all tensorized Quad elements. It should concretely define things common
   * to all Quads. For example, here we define specific methods to attach vertex coordinates to each element, as well
   * as sources and receivers. These last two methods bring up an interesting point however. As you might notice, the
   * mixin approach lets us be perform class nesting to an arbitrary depth, without incurring any runtime cost. So,
   * for example, we could have a super class called something like Element2D, where source and receiver attachments
   * are implemented. However, at one point we need to decide on an optimal level of complexity. Since Quad is the
   * absolute base for all further derivations, this means that attachReceiver() and attachSource() will need to be
   * duplicated for Triangles, and for all other 2D elements. Is this an acceptable level of code duplication, in order
   * to cut down on the class heiarchy depth? Perhaps. Either way, we should discuss this, and also note that following
   * these mixin principles, things are easily modified.
   *
   * To initialize a new Quad element, you'll need at least one template parameter. For instance:
   *
   * \code
   * auto elem = Quad<QuadP1>
   * auto elem = Quad<Acoustic<QuadP1>>
   * auto elem = Quad<Gravity<Attenuation<Acoustic<QuadP1>>>>
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
  const static int mNumDim = 2;
  const static int mNumVtx = 4;

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
  Eigen::Vector2d mElmCtr;

  // Closure mapping.
  Eigen::VectorXi mClsMap;

  // Quadrature parameters.
  Eigen::VectorXd mIntCrdR;
  Eigen::VectorXd mIntCrdS;
  Eigen::VectorXd mIntWgtR;
  Eigen::VectorXd mIntWgtS;

  // Matrix holding gradient information.
  Eigen::MatrixXd mGrd;

  // Material parameters.
  std::map<std::string,Eigen::Vector4d> mPar;

  // Sources and receivers.
  std::vector<std::shared_ptr<Source>> mSrc;
  std::vector<std::shared_ptr<Receiver>> mRec;

 public:

  /** Allocates memory for work arrays, most private variables. */
  Quad<ConcreteShape>(Options options);

  /** Sets up the test function parameters. */
  static Eigen::VectorXd GllPointsForOrder(const int order);
  static Eigen::VectorXi ClosureMappingForOrder(const int order);
  static Eigen::VectorXd GllIntegrationWeightsForOrder(const int order);
  static Eigen::MatrixXd setupGradientOperator(const int order);
  static Eigen::VectorXd interpolateLagrangePolynomials(const double r, const double s, const int order);

  /**
   * Returns an optimized stride along the r direction.
   * @param [in] f Function defined at GLL points.
   * @param [in] s_ind Index of GLL points along s-axis.
   * @param [in] numPtsS Number of integration points along the s-axis.
   * @param [in] numPtsR Number of integration points along the r-axis.
   */
  static Eigen::VectorXd rVectorStride(const Eigen::Ref<const Eigen::VectorXd>& f,
                                       const int s_ind, const int numPtsS,
                                       const int numPtsR);
  /**
   * Returns an optimized stride along the s direction.
   * @param [in] f Function defined at GLL points.
   * @param [in] s_ind Index of GLL points along s-axis.
   * @param [in] numPtsS Number of integration points along the s-axis.
   * @param [in] numPtsR Number of integration points along the r-axis.
   *
   */
  static Eigen::VectorXd sVectorStride(const Eigen::Ref<const Eigen::VectorXd>& f,
                                       const int r_ind, const int numPtsS,
                                       const int numPtsR);

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
   *
   */
  Eigen::VectorXd getDeltaFunctionCoefficients(const double r, const double s);

  /**
   * Given a model, save the material parameters at the element vertices.
   * @param [in] model The model containing the material parameters.
   * @param [in] parameter The parameter to save.
   */
  void attachMaterialProperties(const ExodusModel *model, std::string parameter);

  /**
   * Given some field at the GLL points, interpolate the field to some general point.
   * @param [in] pnt Position in reference coordinates.
   */
  Eigen::MatrixXd interpolateFieldAtPoint(const Eigen::VectorXd &pnt) { return Eigen::MatrixXd(1, 1); }

  // Setters.
  inline void SetNumNew(const PetscInt num) { mElmNum = num; }
  inline void SetVtxCrd(const Eigen::Ref<const Eigen::Matrix<double,4,2>> &v) { mVtxCrd = v; }

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
  std::vector<std::shared_ptr<Source>> Sources() { return mSrc; }

  // Delegates.
  std::tuple<Eigen::VectorXd, Eigen::VectorXd> buildNodalPoints() {
    return ConcreteShape::buildNodalPoints(mIntCrdR, mIntCrdS, mVtxCrd);
  };

};


