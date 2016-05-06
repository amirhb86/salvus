#pragma once

#include <Element.h>
#include <Element/HyperCube/QuadNew.h>
#include <Element/HyperCube/Quad/QuadP1.h>
template <typename T>
class ElementAdapter: public Element, public T {

 public:

  ElementAdapter(Options options): T(options) {};

  virtual std::shared_ptr<Element> clone() const {
    return std::shared_ptr<Element> (new ElementAdapter(*this));
  }

  virtual void attachMaterialProperties(ExodusModel *model) {
    T::attachMaterialProperties(model);
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
  virtual void recordField(const Eigen::MatrixXd &u) {
    T::recordField(u);
  }

  virtual Eigen::MatrixXd computeSourceTerm(double time) {
    return T::computeSourceTerm(time);
  }
  virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement) {
    return T::computeStiffnessTerm(displacement);
  }
  virtual std::vector<std::string> PullElementalFields() const {
    return T::PullElementalFields();
  }
  virtual std::vector<std::string> PushElementalFields() const {
    return T::PushElementalFields();
  }
  virtual Eigen::MatrixXd interpolateFieldAtPoint(const Eigen::VectorXd &pnt) {
    return T::interpolateFieldAtPoint(pnt);
  }
};

template class ElementAdapter<AcousticNew<QuadNew<QuadP1>>>;