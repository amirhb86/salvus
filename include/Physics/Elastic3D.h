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
class Elastic3D: public Shape {
  /**
   * \class Elastic3D
   *
   * \brief Class in charge of handling wave propagation in elastic volumes.
   *
   * This element expects to be templated on "Shape", which refers to a concrete element type.
   * Some examples might be "TensorQuad", or "Generic". Functionality from these derived classes
   * will be required in, for example, the stiffness routine (to calculate the strain).
   */

 private:
  /**** Workspace vectors (allocated in the constructor). ****/
  Eigen::ArrayXd mc11, mc12, mc13, mc22, mc23, mc33, mc44, mc55, mc66, mRho;
//  Eigen::MatrixXd mStiff, mStress, mStrain;

 public:

  /**** Initializers ****/
  Elastic3D<Shape>(std::unique_ptr<Options> const &options);
  std::vector<std::string> PullElementalFields() const;
  std::vector<std::string> PushElementalFields() const;

  /**** Setup functions ****/
  void prepareStiffness();
  void assembleElementMassMatrix(std::unique_ptr<Mesh> const &mesh);
  void attachMaterialProperties(std::unique_ptr<ExodusModel> const &model);
  double CFL_estimate();

  /**** Time loop functions ****/
  Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &u);
  Eigen::MatrixXd computeSourceTerm(const double time);
  Eigen::MatrixXd computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd>& u);
  Eigen::Array<double,Eigen::Dynamic,6> computeStress(const Eigen::Ref<const Eigen::ArrayXd>& strain);
  void recordField(const Eigen::MatrixXd &u) {};

  /**** Test helpers ****/
  void setupEigenfunctionTest(std::unique_ptr<Mesh> const &mesh, std::unique_ptr<Options> const &options) {};
  double checkEigenfunctionTest(std::unique_ptr<Mesh> const &mesh, std::unique_ptr<Options> const &options,
                                const Eigen::Ref<const Eigen::MatrixXd>& u,
                                double time) {return 0;};

};


