#pragma once

#include <ElementNew.h>
#include <Element/HyperCube/QuadNew.h>
#include <Element/HyperCube/Quad/QuadP1.h>
template <typename T>
class ElementAdapter: public ElementNew, public T {

 public:

  ElementAdapter(Options options): T(options) {};

  virtual std::shared_ptr<ElementNew> clone() const {
    return std::shared_ptr<ElementNew> (new ElementAdapter(*this));
  }

  virtual void attachMaterialPropertiesNew(ExodusModel *model) {
    T::attachMaterialPropertiesNew(model);
  }
  virtual void attachVertexCoordinates(DM &distributed_mesh) {
    T::attachVertexCoordinates(distributed_mesh);
  }
  virtual void attachSource(std::vector<std::shared_ptr<Source>> sources) {
    T::attachSource(sources);
  }
  virtual void attachReceiver(std::vector<std::shared_ptr<Receiver>> &receivers) {
    T::attachReceiver(receivers);
  }
  virtual void prepareStiffness() {
    T::prepareStiffness();
  }
  virtual void assembleElementMassMatrix(Mesh *mesh) {
    T::assembleElementMassMatrix(mesh);
  }
  virtual void recordField(const Eigen::Ref<const Eigen::MatrixXd>& field) {
    T::recordField(field);
  }

  virtual Eigen::MatrixXd computeSourceTerm(double time) {
    return T::computeSourceTerm(time);
  }
  virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::Ref<const Eigen::MatrixXd>& u) {
    return T::computeStiffnessTerm(u);
  }
  virtual std::vector<std::string> PullElementalFields() const {
    return T::PullElementalFields();
  }
  virtual std::vector<std::string> PushElementalFields() const {
    return T::PushElementalFields();
  }
  virtual Eigen::MatrixXd interpolateFieldAtPoint(const Eigen::Ref<const Eigen::VectorXd>& pnt) {
    return T::interpolateFieldAtPoint(pnt);
  }
  virtual void setBoundaryConditions(Mesh *mesh) {
    T::setBoundaryConditionsNew(mesh);
  }
  virtual void setupTest(Mesh *mesh, Options options) {
    T::setupEigenfunctionTest(mesh, options);
  }
  virtual double checkTest(Mesh *mesh, Options options, const Eigen::Ref<const Eigen::MatrixXd>& u, double time) {
    return T::checkEigenfunctionTest(mesh, options, u, time);
  }
  virtual void applyDirichletBoundaries(Mesh *mesh, Options options, const std::string &fieldname) {
    return T::applyDirichletBoundaries(mesh, options, fieldname);
  }


  inline void SetNum(const int num) { T::SetNumNew(num); }
  inline int NumDim() const { return T::NumDimNew(); }

  inline int Num() const { return T::ElmNum(); }
  inline int NumIntPnt() const { return T::NumIntPnt(); }
  inline Eigen::VectorXi ClsMap() const { return T::ClosureMap(); }


  inline bool BndElm() const { return T::BndElm(); }
  inline int NumDofVol() const { return T::NumDofVolNew(); }
  inline int NumDofFac() const { return T::NumDofFacNew(); }
  inline int NumDofEdg() const { return T::NumDofEdgNew(); }
  inline int NumDofVtx() const { return T::NumDofVtxNew(); }

};

template class ElementAdapter<AcousticNew<QuadNew<QuadP1>>>;