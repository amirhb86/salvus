#pragma once

#include <Eigen/Dense>
#include <Element/ElementNew.h>

// Derived.
#include <Physics/AcousticNew.h>

extern "C" {
#include <Quad/Autogen/quad_autogen.h>
}

template <typename Derived>
class QuadNew: public Derived {

  /** \class QuadNew
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
   * You can see how the polymorphism is handled by looking at the 'Override' functions below. Basically, at the bottom
   * of all element class hierarchies is a purely static class handling coordinate transforms. This class stores no
   * information itself, it really just provides efficient routines for calculating the Jacobian, interpolation
   * routines, etc. This is nice because it can be composed with anything, but also has the advantage of being purely
   * standalone (i.e. you don't need to initialize it). This makes things like testing very easy. It also means that
   * you need to pass down instance variables (such as the vertex coordinates) when you call these functions from
   * within a class.
   */

private:

  // Static variables.
  const static int mNumDim = 2;
  const static int mNumVtx = 4;

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

  QuadNew<Derived>(Options options);

  static Eigen::VectorXd GllPointsForOrder(const int order);
  static Eigen::VectorXi ClosureMappingForOrder(const int order);
  static Eigen::VectorXd GllIntegrationWeightsForOrder(const int order);
  static Eigen::MatrixXd setupGradientOperator(const int order);

  // Helper functions.
  static Eigen::VectorXd rVectorStride(const Eigen::Ref<const Eigen::VectorXd>& f,
                                       const int s_ind, const int numPtsS,
                                       const int numPtsR) {
    return Eigen::Map<const Eigen::VectorXd> (
        f.data() + s_ind * numPtsS, numPtsR);
  };
  static Eigen::VectorXd sVectorStride(const Eigen::Ref<const Eigen::VectorXd>& f,
                                       const int r_ind, const int numPtsS,
                                       const int numPtsR) {
    return Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>> (
        f.data() + r_ind, numPtsS, Eigen::InnerStride<>(numPtsR));
  };


  void setBoundaryConditionsNew(Mesh *mesh);
  double integrateField(const Eigen::Ref<const Eigen::VectorXd>& field);
  void prepareStiffness();
  void attachVertexCoordinates(DM &distributed_mesh);
  void attachSource(std::vector<std::shared_ptr<Source>> sources);
  void attachReceiver(std::vector<std::shared_ptr<Receiver>> &receivers);
  void applyDirichletBoundaries(Mesh *mesh, Options &options, const std::string &fieldname);

  // Setup.
  void attachMaterialProperties(ExodusModel *model) {};
  void attachMaterialPropertiesNew(const ExodusModel *model, std::string parameter);
  void assembleElementMassMatrix(Mesh *mesh) {};
  void setupTest(Mesh *mesh, Options options);

  double ParAtPnt(const double r, const double s, const std::string &par);

  // Time loop.
  Eigen::MatrixXd computeSourceTerm(double time) { return Eigen::MatrixXd(1, 1); };


  Eigen::MatrixXd interpolateFieldAtPoint(const Eigen::VectorXd &pnt) { return Eigen::MatrixXd(1, 1); }
  void recordField(const Eigen::MatrixXd &u) {};

  // Setters.
  inline void SetNumNew(const PetscInt num) { mElmNum = num; }
  inline void SetVtxCrd(const Eigen::Ref<const Eigen::Matrix<double,4,2>> &v) { mVtxCrd = v; }

  // Getters.
  inline PetscInt ElmNum() const { return mElmNum; }
  inline int NumDimNew() const { return mNumDim; }
  inline int NumIntPnt() const { return mNumIntPnt; }
  inline int NumIntPtsR() { return mNumIntPtsR; }
  inline int NumIntPtsS() { return mNumIntPtsS; }

  inline int NumDofVolNew() const { return mNumDofVol; }
  inline int NumDofFacNew() const { return mNumDofFac; }
  inline int NumDofEdgNew() const { return mNumDofEdg; }
  inline int NumDofVtxNew() const { return mNumDofVtx; }

  inline double IntWgtR(const int ind) { return mIntWgtR[ind]; }
  inline double IntWgtS(const int ind) { return mIntWgtS[ind]; }
  inline double IntCrdR(const int ind) { return mIntCrdR[ind]; }
  inline double IntCrdS(const int ind) { return mIntCrdS[ind]; }
  inline Eigen::VectorXd VecIntWgtR() { return mIntWgtR; }
  inline Eigen::VectorXd VecIntWgtS() { return mIntWgtS; }

  inline bool BndElm() const { return mBndElm; }
  inline Eigen::MatrixXi ClosureMap() const { return mClsMap; }
  inline Eigen::Matrix<double,mNumVtx,mNumDim> VtxCrd() { return mVtxCrd; }
  inline Eigen::VectorXd GrdRow(const int ind) { return mGrd.row(ind); }
  inline Eigen::VectorXd GrdCol(const int ind) { return mGrd.col(ind); }

  // Delegates.
  std::tuple<Eigen::VectorXd, Eigen::VectorXd> buildNodalPoints() {
    return Derived::buildNodalPoints(mIntCrdR, mIntCrdS, mVtxCrd);
  };

};


