#pragma once

#include <Eigen/Dense>
#include <Element/Element.h>

// Physics.
#include <Physics/AcousticNew.h>

using namespace Eigen;

template <typename Physics>
class QuadNew: public Element, Physics {

 private:

  const static int mNumDim = 2;
  const static int mNumVtx = 4;
  Eigen::Matrix<double,mNumVtx,mNumDim> mVtxCrd;

 public:

  std::shared_ptr<Element> clone() const { return std::shared_ptr<Element> (new QuadNew(*this)); }

  void prepareStiffness();

  // Setup.
  void attachMaterialProperties(ExodusModel *model) {};
  void attachVertexCoordinates(DM &distributed_mesh);
  void attachSource(std::vector<std::shared_ptr<Source>> sources) {};
  void attachReceiver(std::vector<std::shared_ptr<Receiver>> &receivers) {};
  void assembleElementMassMatrix(Mesh *mesh) {};

  // Time loop.
  MatrixXd computeSourceTerm(double time) { return MatrixXd(1, 1); };
  MatrixXd computeStiffnessTerm(const MatrixXd &displacement) { return MatrixXd(1, 1); }

  // Utility.
  std::vector<std::string> PullElementalFields() const { return std::vector<std::string> { "Hi" }; }
  std::vector<std::string> PushElementalFields() const { return std::vector<std::string> { "Hi" }; }
  MatrixXd interpolateFieldAtPoint(const VectorXd &pnt) { return MatrixXd(1, 1); }
  void recordField(const MatrixXd &u) {};

};


