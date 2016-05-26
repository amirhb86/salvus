#pragma once

#include <Eigen/Dense>
#include <vector>

class Mesh;
class Options;
class ExodusModel;

template <typename Shape>
class AcousticTri: public Shape {
  /**
   * \class Acoustic2D
   *
   * \brief Class in charge of handling wave propagation in acoustic regions.
   *
   * This element expects to be templated on "Shape", which refers to a concrete element type.
   * Some examples might be "TensorQuad", or "Generic". Functionality from these derived classes
   * will be required in, for example, the stiffness routine (to calculate the strain).
   */

 private:

  /**** Workspace vectors (allocated in the constructor). ****/
  Eigen::VectorXd mVpSquared;
  Eigen::VectorXd mStiff;
  Eigen::VectorXd mSource;
  Eigen::MatrixXd mStress;
  Eigen::MatrixXd mStrain;
  Eigen::MatrixXd mElementStiffnessMatrix;
 public:

  /**** Initializers ****/
  AcousticTri<Shape>(Options options);
  std::vector<std::string> PullElementalFields() const;
  std::vector<std::string> PushElementalFields() const;

  /**** Setup functions ****/
  void prepareStiffness();
  void assembleElementMassMatrix(Mesh *mesh);
  void attachMaterialPropertiesNew(const ExodusModel *model);

  /**** Time loop functions ****/
  Eigen::MatrixXd computeStress(const Eigen::Ref<const Eigen::MatrixXd>& strain);
  Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &u);
  Eigen::MatrixXd computeSourceTerm(const double time);
  void recordField(const Eigen::MatrixXd &u) {};

  /**** Test helpers ****/
  void setupEigenfunctionTest(Mesh *mesh, Options options);
  double checkEigenfunctionTest(Mesh *mesh, Options options,
                                const Eigen::Ref<const Eigen::MatrixXd>& u,
                                double time);

};
