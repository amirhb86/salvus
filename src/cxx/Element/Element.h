#pragma once

#include <Eigen/Dense>
#include <Model/ExodusModel.h>
#include <Source/Source.h>
#include <Mesh/Mesh.h>
#include <mpi.h>
#include <petsc.h>

using namespace Eigen;

class Element {

 protected:

  // Element description.
  int mElmNum;      /** < Element number (on local proc). */
  int mPlyOrd;      /** < Polynomial order. */
  int mNumDim;      /** < Num of dimensions. */
  int mNumDofVtx;   /** < Num DOF per vertex. */
  int mNumDofEdg;   /** < Num DOF per edge. */
  int mNumDofFac;   /** < Num DOF per face. */
  int mNumDofVol;   /** < Num DOF per volume. */
  int mNumIntPnt;   /** < Total number of integration points. */

  // Element properties.
  VectorXi mClsMap; /** < Mapping from our element closure numbering to PETSc's. */
  VectorXd mMssMat; /** < Elemental mass matrix. */
  MatrixXd mVtxCrd; /** < Vertex coordinates ordered as above. r0->x,r1->y,r2->z. */
  MatrixXd mElmCtr; /** < Coordinates of the element center. r0->x,r1->y,r2->z. */

  // Boundary information.
  bool mBndElm;     /** < Is element on a boundary. */
  std::map<std::string,
           std::vector<int>> mBnd;  /** < Map relating the type of boundary to the Petsc edge number */

  // Other objects.
  std::vector<Source *> mSrc; /** < Vector of sources belonging to the element. */

 public:

  /**
   * Factory return the proper element physics based on the command line options.
   * @return Some derived element class.
   */
  static Element *factory(Options options);

  /**
   * Copy constructor.
   * Returns a copy. Use this once the reference element is set up via the constructor, to allocate space for
   * all the unique elements on a processor.
   */
  virtual Element *clone() const = 0;

  /**
   * Utility function to integrate a field over the element. This could probably be made static, but for now I'm
   * just using it to check some values.
   * @param [in] field The field which to integrate, defined on each of the gll points.
   * @returns The scalar value of the field integrated over the element.
   */
  virtual double integrateField(const Eigen::VectorXd &field) = 0;

  /**
   * Queries the passed DM for the vertex coordinates of the specific element. These coordinates are saved
   * in mVertexCoordiantes.
   * @param [in] distributed_mesh PETSc DM object.
   *
   */
  virtual void attachVertexCoordinates(DM &distributed_mesh) = 0;

  /**
   * Attach source.
   * Given a vector of abstract source objects, this function will query each for its spatial location. After
   * performing a convex hull test, it will perform a quick inverse problem to determine the position of any sources
   * within each element in reference coordinates. These reference coordinates are then saved in the source object.
   * References to any sources which lie within the element are saved in the mSrc vector.
   * @param [in] sources A vector of all the sources defined for a simulation run.
   */
  virtual void attachSource(std::vector<Source *> sources) = 0;

  /**
   * Build the elemental stiffness matrix.
   * If necessary, construct K on each element. This requires quantities like the material
   * parameters to already be defined. As well, in the case of the tensorized basis, we take
   * advantage of the sparsity structure of the stiffness matrix to skip this step altogether.
   */
  virtual void prepareStiffness() = 0;

  /**
   * TODO: Get rid of mesh as a parameter.
   * Build the elemental mass matrix.
   * Populate the coefficients of the mass matrix. When this function has returned, the mass matrix
   * should be ready to scatter to the global degrees of freedom.
   */
  virtual void assembleElementMassMatrix(Mesh *mesh) = 0;

  /**
   * TODO: Make this general over dimension.
   * Builds nodal coordinates (x,(y),z) on all mesh degrees of freedom.
   * @param mesh [in] The mesh.
   */
  virtual std::tuple<Eigen::VectorXd, Eigen::VectorXd> buildNodalPoints() = 0;

  /**
   * TODO: A good way to represent a completely general source.
   * Compute the right hand side.
   * Given a certain time, compute the value of the external force.
   */
  virtual MatrixXd computeSourceTerm(double time) = 0;

  /**
   * Interpolate material parameters to element edge.
   * A lightweight (in memory) way to carry around material parameters is to just store them
   * at the edges of elements. On the fly, during the time loop, interpolation to integration points
   * can be performed by a quick dot product.
   * @param model [in] A class which returns a model model parameter at a queried point.
   */
  virtual void interpolateMaterialProperties(ExodusModel *model) = 0;

  // Attribute sets.
  void SetNum(int element_number) { mElmNum = element_number; }
  void SetVtxCrd(MatrixXd coord) { mVtxCrd = coord; }

  // Attribute gets.
  inline int Num() const { return mElmNum; }
  inline int NumDim() const { return mNumDim; }
  inline int NumDofEdg() const { return mNumDofEdg; }
  inline int NumDofFac() const { return mNumDofFac; }
  inline int NumDofVtx() const { return mNumDofVtx; }
  inline int NumIntPnt() const { return mNumIntPnt; }
  inline bool BndElm() const { return mBndElm; }
  VectorXi ClsMap() { return mClsMap; }
  MatrixXd VtxCrd() { return mVtxCrd; }

  /**
   * Compute the action of the stiffness matrix on a field.
   *
   */
  virtual MatrixXd computeStiffnessTerm(const MatrixXd &displacement) = 0;

  /**
   * Figure out which dofs (if any) are on the boundary.
   */
  virtual void setBoundaryConditions(Mesh *mesh);

  /**
   * Apply boundary condition (dirichlet, absorbing, fluid-elastic coupling, etc)
   */
  virtual void applyBoundaryConditions(Mesh *mesh,
                                       Options &options,
                                       std::string fieldname);


  /**
   * Returns the fields which are required on the element level.
   */
  virtual std::vector<std::string> PullElementalFields() const = 0;

  /**
   * Returns the fields which the element assembles into the global dofs.
   */
  virtual std::vector<std::string> PushElementalFields() const = 0;

  // Testing, which is usually implemented via an initial condition
  virtual void setupTest(Mesh *mesh, Options options);
  virtual double checkTest(Mesh *mesh, Options options, const Eigen::MatrixXd &displacement, double time);

};
