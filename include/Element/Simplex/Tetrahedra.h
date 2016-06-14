#pragma once

// stl.
#include <map>
#include <vector>
#include <iostream>

// 3rd party.
#include <Eigen/Dense>

// forward decl.
class Mesh;
class Source;
class Options;
class Receiver;
class ExodusModel;

template <typename ConcreteShape>
class Tetrahedra: public ConcreteShape {
  /**
   * Base class of an abstract four node tetrahedra. The reference element is set up as below.
   * 
   *              (v3)
   *                +--
   *                | \----
   *                |   \--\----
   *                |      \--  \----
   *                |         \-     \---  (v1)
   *                |           \--      \--
   *  ^             |              \--/---  \
   *  |             |             /---\-     \
   *  | (t)         |         /---      \--   \
   *  |             |     /---             \-- \
   *  |      (s)    | /---                    \-\
   *  |       /-    +---------------------------\--
   *  |    /--     (v0)                        (v2)
   *  | /--
   *  +-----------------> (r)
   * 
   * Faces are organized like: bottom (r-s), left (s-t), front (r-t), right (r-s-t)
   * More precisely, faces are composed of vertices and ordered
   * 0: 0 1 2
   * 1: 0 3 1
   * 2: 0 2 3
   * 3: 2 1 3
   * Edge ordering
   * (v0,v1)
   * (v1,v2)
   * (v2,v0)
   * (v0,v3)
   * (v3,v1)
   * (v2,v3)
   */
    
 private:

  /*****************************************************************************
   * STATIC MEMBERS. THESE VARIABLES AND FUNCTIONS SHOULD APPLY TO ALL ELEMENTS.
   *****************************************************************************/

  // Static variables.
  const static int mNumDim = 3;
  const static int mNumVtx = 4;

  // Instance variables.
  PetscInt mElmNum;
  int mPlyOrd;
  int mNumIntPnt;
  int mNumDofVtx;
  int mNumDofEdg;
  int mNumDofFac;
  int mNumDofVol;

  static Eigen::MatrixXd mGradientOperator;
  /** < Derivative of shape function n (col) at pos. m (row) */
  static Eigen::VectorXd mIntegrationWeights;
  /** < Integration weights along epsilon direction. */
  static Eigen::VectorXd mIntegrationCoordinates_r;
  /** < Nodal location directions */
  static Eigen::VectorXd mIntegrationCoordinates_s;  
  static Eigen::VectorXd mIntegrationCoordinates_t;  

  /** r,s,t-derivative of Lagrange poly Phi [row: phi_i, col: phi @ nth point]*/
  static Eigen::MatrixXd mGradientPhi_dr;
  static Eigen::MatrixXd mGradientPhi_ds;
  static Eigen::MatrixXd mGradientPhi_dt;
  static Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> mGradientPhi_dr_t;
  static Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> mGradientPhi_ds_t;
  static Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> mGradientPhi_dt_t;
  
  /**********************************************************************************
   * OBJECT MEMBERS. THESE VARIABLES AND FUNCTIONS SHOULD APPLY TO SPECIFIC ELEMENTS.
   ***********************************************************************************/

  Eigen::MatrixXd mGradientPhi_dx;
  Eigen::MatrixXd mGradientPhi_dy;
  Eigen::MatrixXd mGradientPhi_dz;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> mGradientPhi_dx_t;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> mGradientPhi_dy_t;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> mGradientPhi_dz_t;

  Eigen::MatrixXd mWiDPhi_x;
  Eigen::MatrixXd mWiDPhi_y;
  Eigen::MatrixXd mWiDPhi_z;
  
  // Vertex coordinates.
  Eigen::Matrix<double,mNumVtx,mNumDim> mVtxCrd;
  
  // Element center.
  Eigen::Vector3d mElmCtr;

  // Closure mapping.
  Eigen::VectorXi mClsMap;
  
  // Workspace.
  double mDetJac;
  Eigen::VectorXd mParWork;
  Eigen::VectorXd mStiffWork;
  Eigen::MatrixXd mGradWork;
  Eigen::VectorXd mGrad_r;
  Eigen::VectorXd mGrad_s;
  Eigen::VectorXd mGrad_t;
  
  // On Boundary.
  bool mBndElm;
  std::map<std::string,std::vector<int>> mBnd;

  // Material parameters.
  std::map<std::string,Eigen::Vector4d> mPar;

  // Sources and receivers.
  std::vector<std::shared_ptr<Source>> mSrc;
  std::vector<std::shared_ptr<Receiver>> mRec;
  
  // precomputed element stiffness matrix (with velocities)
  Eigen::MatrixXd mElementStiffnessMatrix;

  Eigen::Matrix3d mInvJac;
  Eigen::Matrix3d mInvJacT;
  Eigen::Matrix3d mInvJacT_x_invJac;
  
 public:

  /**
   * Factory return the proper element physics based on the command line options.
   * @return Some derived element class.
   */
  static Tetrahedra *factory(Options options);

  /**
   * Constructor.
   * Sets quantities such as number of dofs, among other things, from the options class.
   * @param [in] options Populated options class.
   */
  Tetrahedra(Options options);

  /**
   * Returns the gll locations for a given polynomial order.
   * @param [in] order The polynmomial order.
   * @returns tuple of quadrature points (r,s).
   */
  static std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> QuadraturePoints(const int order);

  /**
   * Returns the quadrature intergration weights for a polynomial order.
   * @param [in] order The polynomial order.
   * @returns Vector of quadrature weights.
   * TODO: Move to autogenerated code.
   */
  static Eigen::VectorXd QuadratureIntegrationWeights(const int order);

  /**
   * Returns the mapping from the PETSc to Salvus closure.
   * @param [in] order The polynomial order.
   * @param [in] dimension Element dimension.
   * @returns Vector containing the closure mapping (field(closure(i)) = petscField(i))
   */
  Eigen::VectorXi ClosureMapping(const int order, const int dimension,
                                 DM &distributed_mesh);

  static Eigen::Vector4d interpolateAtPoint(double r, double s, double t);
  Eigen::MatrixXd interpolateFieldAtPoint(const Eigen::VectorXd &pnt) {
    std::cerr << "ERROR: interpolateFieldAtPoint not implemented\n"; exit(1);
    return Eigen::Matrix<double, -1, -1, 0, -1, -1>();
  }

  void recordField(const Eigen::MatrixXd &u) { std::cerr << "ERROR: recordField not implemented\n"; exit(1); };

  // currently not possible
  // /**
  //  * Returns the face mapping from the PETSc to Salvus closure.
  //  * @param [in] order The polynomial order.
  //  * @param [in] dimension Element dimension.
  //  * @returns Vector containing the closure mapping (field(closure(i)) = petscField(i))
  //  */
  // static Eigen::VectorXi FaceClosureMapping(const int order, const int dimension);

  // Attribute gets.

  // virtual Eigen::VectorXi GetFaceClosureMapping() { return mFaceClosureMapping; }

  // virtual Eigen::MatrixXd VtxCrd() { return mVtxCrd; }

  /********************************************************
   ********* Inherited methods from Element 2D ************
   ********************************************************/

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
   * Attaches a material parameter to the vertices on the current element.
   * Given an exodus model object, we use a kD-tree to find the closest parameter to a vertex. In practice, this
   * closest parameter is often exactly coincident with the vertex, as we use the same exodus model for our mesh
   * as we do for our parameters.
   * @param [in] model An exodus model object.
   * @param [in] parameter_name The name of the field to be added (i.e. velocity, c11).
   */
  void attachMaterialProperties(const ExodusModel *model,
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
   *
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
   * Atttach receiver.
   * Given a vector of abstract receiver objects, this function will query each for its spatial location. After
   * performing a convex hull test, it will perform a quick inverse problem to determine the position of any
   * sources within each element in reference coordiantes. These reference coordinates are then saved in the
   * receiver object. References to any receivers which lie within the element are saved in the mRec vector.
   * @param [in] receivers A vector of all the receivers defined for a simulation run.
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
  Eigen::VectorXd getDeltaFunctionCoefficients(double r, double s, double t) { std::cerr << "ERROR: Not implemented"; exit(1); }

  Eigen::VectorXd computeStiffnessFull(const Eigen::Ref<const Eigen::VectorXd>&  field,
                                       const Eigen::Ref<const Eigen::VectorXd>& vp2);
  
  Eigen::MatrixXd buildStiffnessMatrix(const Eigen::Ref<const Eigen::VectorXd>& vp2);

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

  double CFL_constant();
  
  /** Return the estimated element radius
   * @return The CFL estimate
   */
  double estimatedElementRadius();
  
  /**
   * Simple function to set the (remembered) element number.
   */
  void SetNum(int element_number) { mElmNum = element_number; }
  inline void SetNumNew(const PetscInt num) { mElmNum = num; }
  inline void SetVtxCrd(const Eigen::Ref<const Eigen::Matrix<double,mNumVtx,mNumDim>> &v) { mVtxCrd = v; }

  inline PetscInt ElmNum() const { return mElmNum; }
  inline bool BndElm() const { return mBndElm; }
  inline int NumDim() const { return mNumDim; }
  inline int NumIntPnt() const { return mNumIntPnt; }
  inline int NumDofVol() const { return mNumDofVol; }
  inline int NumDofFac() const { return mNumDofFac; }
  inline int NumDofEdg() const { return mNumDofEdg; }
  inline int NumDofVtx() const { return mNumDofVtx; }
  inline Eigen::MatrixXi ClsMap() const { return mClsMap; }
  inline Eigen::MatrixXd VtxCrd() const { return mVtxCrd; }
  std::vector<std::shared_ptr<Source>> Sources() { return mSrc; }

  /**
   * Builds nodal coordinates (x,z) on all mesh degrees of freedom.
   * @param mesh [in] The mesh.
   */
  // Delegates.
  std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> buildNodalPoints() {
    return ConcreteShape::buildNodalPoints(mIntegrationCoordinates_r,
                                           mIntegrationCoordinates_s,
                                           mIntegrationCoordinates_t,
                                           mVtxCrd);
  };

  
};


