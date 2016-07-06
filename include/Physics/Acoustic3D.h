#pragma once

// stl.
#include <vector>
#include <memory>

// 3rd party.
#include <Eigen/Dense>

// forward decl.
class Mesh;
class Options;
class ExodusModel;

template <typename Shape>
class Acoustic3D: public Shape {
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

 public:

  /**** Initializers ****/
  Acoustic3D<Shape>(std::unique_ptr<Options> const &options);
  std::vector<std::string> PullElementalFields() const;
  std::vector<std::string> PushElementalFields() const;

  /**** Setup functions ****/
  void prepareStiffness();
  Eigen::MatrixXd assembleElementMassMatrix();
  void attachMaterialProperties(std::unique_ptr<ExodusModel> const &model);
  double CFL_estimate();
  
  
  /**** Time loop functions ****/
  Eigen::MatrixXd computeStress(const Eigen::Ref<const Eigen::MatrixXd>& strain);
  Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &u);
  Eigen::MatrixXd computeSourceTerm(const double time);
  Eigen::MatrixXd computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd>& u);
  void recordField(const Eigen::MatrixXd &u) {};

  /**** Test helpers ****/
  void setupEigenfunctionTest(std::unique_ptr<Mesh> const &mesh, std::unique_ptr<Options> const &options);
  double checkEigenfunctionTest(std::unique_ptr<Mesh> const &mesh, std::unique_ptr<Options> const &options,
                                const Eigen::Ref<const Eigen::MatrixXd>& u,
                                double time);

};

