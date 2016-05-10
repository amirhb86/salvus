#pragma once

#include <Eigen/Dense>
#include <Element/Element.h>
#include <vector>

template <typename Shape>
class AcousticNew: public Shape {

 private:

  Eigen::VectorXd mStiff;
  Eigen::VectorXd mDetJac;
  Eigen::VectorXd mVelGrd;
  Eigen::MatrixXd mStress;
  Eigen::MatrixXd mStrain;

 public:

  // Delegates.
  AcousticNew<Shape>(Options options);

  void attachMaterialPropertiesNew(const ExodusModel *model);
  std::vector<std::string> PullElementalFields() const;
  std::vector<std::string> PushElementalFields() const;

  void setupEigenfunctionTest(Mesh *mesh, Options options);
  double checkEigenfunctionTest(Mesh *mesh, Options options,
                                const Eigen::Ref<const Eigen::MatrixXd>& u, double time);
  void assembleElementMassMatrix(Mesh *mesh);
  Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement);

};

