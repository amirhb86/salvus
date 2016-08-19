#pragma once

// stl.
#include <vector>

// 3rd party.
#include <Eigen/Dense>
#include <Utilities/Types.h>

// forward decl.
class Mesh;
class Options;
class ExodusModel;

template <typename Shape>
class Scalar: public Shape {
  /**
   * \class Scalar
   *
   * \brief Class in charge of handling wave propagation in acoustic regions.
   *
   * This element expects to be templated on "Shape", which refers to a concrete element type.
   * Some examples might be "TensorQuad", or "Generic". Functionality from these derived classes
   * will be required in, for example, the stiffness routine (to calculate the strain).
   */

 private:

  /**** Workspace vectors (allocated in the constructor). ****/
  RealVec mVpSquared;
  RealVec mStiff;
  RealVec mSource;
  RealMat mStress;
  RealMat mStrain;

 public:

  /**** Initializers ****/
  Scalar<Shape>(std::unique_ptr<Options> const &options);
  ~Scalar<Shape>() {};
  std::vector<std::string> PullElementalFields() const;
  std::vector<std::string> PushElementalFields() const;

  /**** Setup functions ****/  
  RealMat assembleElementMassMatrix();
  void attachMaterialProperties(std::unique_ptr<ExodusModel> const &model);

  /**** Time loop functions ****/
  RealMat computeStress(const Eigen::Ref<const RealMat>& strain);
  RealMat computeStiffnessTerm(const Eigen::Ref<const RealMat>& u);
  RealMat computeSurfaceIntegral(const Eigen::Ref<const RealMat>& u);
  RealMat computeSourceTerm(const double time);
  void recordField(const Eigen::Ref<const RealMat>& u) {};

  const static std::string Name() { return "Scalar_" + Shape::Name(); }

};

