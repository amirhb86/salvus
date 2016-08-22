#pragma once

// stl.
#include <vector>

// 3rd party.
#include <Eigen/Dense>

// forward decl.
class Mesh;
class Options;
class ExodusModel;

template <typename Shape>
class Elastic2D: public Shape {
  /**
   * \class Elastic2D
   *
   * \brief Class in charge of handling wave propagation in elastic regions.
   *
   * This element expects to be templated on "Shape", which refers to a concrete element type.
   * Some examples might be "TensorQuad", or "Generic". Functionality from these derived classes
   * will be required in, for example, the stiffness routine (to calculate the strain).
   */

 private:

  /**** Workspace vectors (allocated in the constructor). ****/
  Eigen::VectorXd mc11, mc12, mc13, mc22, mc23, mc33;
  Eigen::MatrixXd mStiff;
  Eigen::MatrixXd mStress;
  Eigen::MatrixXd mStrain;

 public:

  /**** Initializers ****/
  Elastic2D<Shape>(std::unique_ptr<Options> const &options);
  std::vector<std::string> PullElementalFields() const;
  std::vector<std::string> PushElementalFields() const;

  /**** Setup functions ****/
  void prepareStiffness() {};
  Eigen::MatrixXd assembleElementMassMatrix();
  void attachMaterialProperties(std::unique_ptr<ExodusModel> const &model);
  double CFL_estimate() { return 1; }
  
  /**** Time loop functions ****/
  Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &u);
  Eigen::MatrixXd computeSourceTerm(const double time, const PetscInt time_idx);
  Eigen::MatrixXd computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd>& u);
  Eigen::MatrixXd computeStress(const Eigen::Ref<const Eigen::MatrixXd>& strain);
  void recordField(const Eigen::MatrixXd &u) {};

  /**** Test helpers ****/
  const static std::string Name() { return "Elastic2D_" + Shape::Name(); }

};

