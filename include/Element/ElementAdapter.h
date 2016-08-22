#pragma once

#include <Element/Element.h>
#include <Utilities/Options.h>

template <typename T>
class ElementAdapter: public Element, public T {
  /**
   * \class ElementAdapter
   *
   * \brief Class which serves to delegate the abstract Element class to concrete implementations.
   *
   * This class serves as a bridge between the abstract interface defined in Element, and the generic
   * mixin templates defined by T. You see below that nothing concrete really happens here, we just take
   * this opportunity to delegate functionality to the templates. In practice, these delegation functions
   * function similar to virtual functions (that's a mouthful!). Except of course we're not storing pointers
   * to some function, but instead generating code at compile time by combination with the template
   * parameter T.
   *
   * At the end of this file, we define all the concrete types we want. For example
   *
   * \code
   * typedef class ElementAdapter<AcousticNew<QuadNew<QuadP1>>> AcousticQuadP1
   * \endcode
   *
   * defines a new type AcousticQuadP1 that contains all the properly delegated functionality needed
   * to a finite element of that type. Since this type also inherits from the interface class Element,
   * it can be treated as any other generic element, and be (for example) used in a element-type stl
   * container.
   */

 public:

  /**
   * Delegated constructor.
   * Here we (must) implicitly call the templated Base-class constructor. Since the current
   * class is really still just an interface - adapter, there is nothing to initialize here.
   */
  ElementAdapter(std::unique_ptr<Options> const &options): T(options) { };

  /** @name Management.
   * These methods are mainly responsible for memory management, i.e. the creation and desctruction
   * of individual elements.
   */

  /** @name Element setup.
   * These methods are responsible for setting up each individual element for use within a time loop.
   * Basically, we are building up a complex function object.
   */
  ///@{
  /** Construct the mass matrix on the element, and sum into global DOF.
   * @param [in/out] mesh Mesh instance. Global and local vectors modified.
   */
  virtual Eigen::MatrixXd assembleElementMassMatrix() {
     return T::assembleElementMassMatrix();
  }
  /** Attach material parameters to the element given some model.
   * @param [in] model Model instance.
   */
  virtual void attachMaterialProperties(std::unique_ptr<ExodusModel> const &model) {
    T::attachMaterialProperties(model);
  }
  /** Attach receivers to the element (if required).
   * @param [in/out] receivers Vector of all receivers in the model. Receiver reference coordinates are attached.
   */
  virtual bool attachReceiver(std::unique_ptr<Receiver> &receiver, const bool finalize) {
    return T::attachReceiver(receiver, finalize);
  }
  /** Attach sources to the element (if required).
   * @param [in] sources Vector of all sources in the model.
   */
  virtual bool attachSource(std::unique_ptr<Source> &source, const bool finalize) {
    return T::attachSource(source, finalize);
  }
  /** Attach vertex coordinates to the element.
   * @param [in] distributed_mesh The parallel DM provided by PETSc.
   */
  virtual void attachVertexCoordinates(std::unique_ptr<Mesh> const &mesh) {
    T::attachVertexCoordinates(mesh);
  }

  /** Precompute any terms needed on the element level, e.g.,
      jacobians, velocities at nodes, stiffness matrices for triangles
      and tetrahedra.
  */
  virtual void precomputeElementTerms() {
    T::precomputeElementTerms();
  }
  ///@}

  /** @name Time loop (pure functions).
   * These functions are called from within the time loop. At this point, the element is a fully constructed
   * function.
   */
  ///@{
  /** Returns the interpolated source for a given time.
   * @ param [in] time Simulation time.
   * @ param [in] time_idx Simulation time index.
   */
  virtual Eigen::MatrixXd computeSourceTerm(const double time, const PetscInt time_idx) {
    return T::computeSourceTerm(time, time_idx);
  }
  /** Returns the action of the stiffness matrix applied to some field (probably displacement).
   * @param [in] u Displacement field.
   */
  virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::Ref<const Eigen::MatrixXd>& u) {
    return T::computeStiffnessTerm(u);
  }
  /** Computes the surface integral over an element. Note that this is usually zero. */
  virtual Eigen::MatrixXd computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd>& u) {
    return T::computeSurfaceIntegral(u);
  };
  /** Returns the fields which are required from the global DOFs for local operation */
  virtual std::vector<std::string> PullElementalFields() const {
    return T::PullElementalFields();
  }
  /** Returns the fields from the global DOFs into which we will sum */
  virtual std::vector<std::string> PushElementalFields() const {
    return T::PushElementalFields();
  }
  /* TODO: Check if the following function is in the right place. */
  /** Returns the (real-space) lagrange polynomials evaluated at some point. */
  virtual Eigen::MatrixXd interpolateFieldAtPoint(const Eigen::Ref<const Eigen::VectorXd>& pnt) {
    return T::interpolateFieldAtPoint(pnt);
  }
  ///@}

  /** @name Time loop (functions with side effects).
   * These functions are also called from within the time loop, but they have a roll which is slightly more
   * interesting than just transforming input into output. Perhaps we may find a way to get rid of them.
   */
  ///@{
  /** Sets the mesh boundaries to a certain value (i.e. Dirichlet, Neumann, ...). Note that this only refers to
   * the actualy physical boundaries of the mesh. Internal material boundaries are handled implicitly by the
   * elements.
   * @param [in/out] mesh Mesh instance. Local/global boundary values are modified.
   */
  virtual void setBoundaryConditions(std::unique_ptr<Mesh> const &mesh) {
    T::setBoundaryConditions(mesh);
  }
  /** Triggers a record on all receivers which belong to a certain element.
   * @param [in] field The field to record.
   */
  virtual void recordField(const Eigen::Ref<const Eigen::MatrixXd>& field) {
    T::recordField(field);
  }

  /** @name Setters/Getters.
   * The time loop occasionally needs access to some instance variables stored in the derived classes. These functions
   * are meant to return directly from those derived classes, so that no instance variables are stored in the interface
   * itself.
   */
  ///@{
  /** Set the element number of the local partition. */
  inline void SetNum(const int num) { T::SetNumNew(num); }
  /** Is this element on a mesh boundary. */
  inline bool BndElm() const { return T::BndElm(); }
  /** What is this elements number on the local partition. */
  inline int Num() const { return T::ElmNum(); }
  /** What is the dimension of this element. */
  inline int NumDim() const { return T::NumDim(); }
  /** How many degrees of freedom in this element's volume. */
  inline int NumDofVol() const { return T::NumDofVol(); }
  /** How many degrees of freedom in this element's face. */
  inline int NumDofFac() const { return T::NumDofFac(); }
  /** How many degrees of freedom on this element's edges. */
  inline int NumDofEdg() const { return T::NumDofEdg(); }
  /** How many degrees of freedom on this element's verticies. */
  inline int NumDofVtx() const { return T::NumDofVtx(); }
  /** How many integration point on this element. */
  inline int NumIntPnt() const { return T::NumIntPnt(); }
  /** Element vertices. */
  inline Eigen::MatrixXd VtxCrd() const { return T::VtxCrd(); }
  /** What is the closure map on this element. */
  virtual inline Eigen::VectorXi ClsMap() const { return T::ClsMap(); }
  /** Vertex coordinates of this element. */
  virtual inline int PlyOrd() const { return T::PlyOrd(); }
  /** What type of element am I? */
  inline std::string Name() const  { return T::Name(); }
  ///@}
};
