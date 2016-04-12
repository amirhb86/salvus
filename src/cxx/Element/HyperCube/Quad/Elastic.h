#pragma once
#include <Eigen/Dense>
#include "../Quad.h"

class Elastic: public Quad {

  static const std::vector<std::string> mElementalFields;
  Eigen::Vector4d mC11AtVertices;
  Eigen::Vector4d mC13AtVertices;
  Eigen::Vector4d mC15AtVertices;
  Eigen::Vector4d mC31AtVertices;
  Eigen::Vector4d mC33AtVertices;
  Eigen::Vector4d mC35AtVertices;
  Eigen::Vector4d mC51AtVertices;
  Eigen::Vector4d mC53AtVertices;
  Eigen::Vector4d mC55AtVertices;
  Eigen::Vector4d mRhoAtVertices;

  Eigen::Matrix<double, 2, Eigen::Dynamic> mElementStrain;

 public:

  Elastic(Options options);
  ~Elastic() {};
  std::shared_ptr<Element> clone() const { return std::shared_ptr<Element> (new Elastic(*this)); }

  virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement);
  virtual Eigen::MatrixXd computeSourceTerm(double time);
  virtual void computeSurfaceTerm();
  virtual void assembleElementMassMatrix(Mesh *mesh);
  virtual void attachMaterialProperties(ExodusModel *model);
  virtual std::vector<std::string> PullElementalFields() const { return {"ux", "uy"}; }
  virtual std::vector<std::string> PushElementalFields() const { return {"ax", "ay"}; }
  void prepareStiffness() { /* nothing to be done for quad */ }

};

