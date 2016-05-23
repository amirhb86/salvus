#pragma once

#include <Eigen/Dense>
#include <vector>

class Mesh;
class Options;
class ExodusModel;

template <typename Element>
class Elastic3D {
  /**
   * \class Elastic3D.
   *
   * \brief Class responsible for 3D fully anisotropic wave propagation.
   *
   * This element class handles the impementation details of the 3D elastic wave equation on generic element
   * types. We formulate our equation in Voigt notation, taking advantage of the symmetric character of the elastic
   * tensor.
   */

 private:

  /**** Workspace vectors (allocated in the constructor). ****/
  Eigen::VectorXd mc11, mc12, mc13, mc14, mc15, mc16;
  Eigen::VectorXd mc22, mc23, mc24, mc25, mc26;
  Eigen::VectorXd mc33, mc34, mc35, mc36;
  Eigen::VectorXd mc44, mc45, mc46;
  Eigen::VectorXd mc55, mc56;
  Eigen::VectorXd mc66;
  Eigen::MatrixXd mStiff, mStress, mStrain;

 public:

  /**** Initializers ****/
  Elastic3D<Element>(Options options);
  std::vector<std::string> PullElementalFields() const;
  std::vector<std::string> PushElementalFields() const;

  /**** Setup functions ****/
  void prepareStiffness() {};
  void assembleElementMassMatrix(Mesh *mesh);
  void attachMaterialPropertiesNew(const ExodusModel *model);

  /**** Time loop functions ****/
  Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &u);
  Eigen::MatrixXd computeSourceTerm(const double time);
  void recordField(const Eigen::MatrixXd &u) {};



};


