#include <Mesh/Mesh.h>
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>
#include <Physics/Elastic2D.h>
#include <Source/Source.h>
#include <Utilities/Types.h>

using namespace Eigen;

template <typename Element>
Elastic2D<Element>::Elastic2D(std::unique_ptr<Options> const &options): Element(options) {

  int num_grad_cmps = 3;

  // Allocate all work arrays.
  mc11.setZero(Element::NumIntPnt());
  mc12.setZero(Element::NumIntPnt());
  mc13.setZero(Element::NumIntPnt());
  mc22.setZero(Element::NumIntPnt());
  mc23.setZero(Element::NumIntPnt());
  mc33.setZero(Element::NumIntPnt());
  mStiff.setZero(Element::NumIntPnt(), Element::NumDim());
  mStress.setZero(Element::NumIntPnt(), num_grad_cmps);
  mStrain.setZero(Element::NumIntPnt(), num_grad_cmps+1);

}

template <typename Element>
void Elastic2D<Element>::attachMaterialProperties(std::unique_ptr<ExodusModel> const &model) {
  Element::attachMaterialProperties(model, "RHO");
  Element::attachMaterialProperties(model, "C11");
//  Element::attachMaterialProperties(model, "C12");
  Element::attachMaterialProperties(model, "C13");
  Element::attachMaterialProperties(model, "C33");
//  Element::attachMaterialProperties(model, "C23");
  Element::attachMaterialProperties(model, "C55");
}

template <typename Element>
std::vector<std::string> Elastic2D<Element>::PullElementalFields() const { return { "ux", "uy"}; }

template <typename Element>
std::vector<std::string> Elastic2D<Element>::PushElementalFields() const { return { "ax", "ay"}; }

template <typename Element>
MatrixXd Elastic2D<Element>::assembleElementMassMatrix() {

  return Element::applyTestAndIntegrate(Element::ParAtIntPts("RHO"));

}

template <typename Element>
MatrixXd Elastic2D<Element>::computeStress(const Eigen::Ref<const Eigen::MatrixXd> &strain) {

  /* TODO: PROPER CONVETION!!!! */
  /* mc13 & mc23 are currently zero */
  mc11 = Element::ParAtIntPts("C11");
  mc12 = Element::ParAtIntPts("C13");
  mc22 = Element::ParAtIntPts("C33");
  mc33 = Element::ParAtIntPts("C55");

  Matrix<double,Dynamic,3> stress(Element::NumIntPnt(), 3);
  VectorXd uxy_plus_uyx = strain.col(1) + strain.col(2);

  stress.col(0) =
    mc11.array().cwiseProduct(strain.col(0).array()) +
    mc12.array().cwiseProduct(strain.col(3).array()) +
    mc13.array().cwiseProduct(uxy_plus_uyx.array());

  stress.col(1) =
      mc12.array().cwiseProduct(strain.col(0).array()) +
      mc22.array().cwiseProduct(strain.col(3).array()) +
      mc23.array().cwiseProduct(uxy_plus_uyx.array());

  stress.col(2) =
      mc13.array().cwiseProduct(strain.col(0).array()) +
      mc23.array().cwiseProduct(strain.col(3).array()) +
      mc33.array().cwiseProduct(uxy_plus_uyx.array());

  return stress;

}

template <typename Element>
MatrixXd Elastic2D<Element>::computeStiffnessTerm(const Eigen::MatrixXd &u) {

  // strain ux_x, ux_y, uy_x, uy_y.
  EIGEN_ASM_COMMENT("BEGIN_GRADIENT");
  mStrain.leftCols<2>()  = Element::computeGradient(u.col(0));
  EIGEN_ASM_COMMENT("END_GRADIENT");
  mStrain.rightCols<2>() = Element::computeGradient(u.col(1));

  // compute stress from strain.
  mStress = computeStress(mStrain);

  // temporary matrix to hold directional stresses.
  Matrix<double,Dynamic,2> temp_stress(Element::NumIntPnt(), 2);

  // compute stiffness.
  temp_stress.col(0) = mStress.col(0); temp_stress.col(1) = mStress.col(2);
  mStiff.col(0) = Element::applyGradTestAndIntegrate(temp_stress);
  temp_stress.col(0) = mStress.col(2); temp_stress.col(1) = mStress.col(1);
  mStiff.col(1) = Element::applyGradTestAndIntegrate(temp_stress);

  return mStiff;

}

template <typename Element>
MatrixXd Elastic2D<Element>::computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd> &u) {
  return Eigen::MatrixXd::Zero(Element::NumIntPnt(), Element::NumDim());
}

template <typename Element>
MatrixXd Elastic2D<Element>::computeSourceTerm(const double time, const PetscInt time_idx) {
  RealMat s = RealMat::Zero(Element::NumIntPnt(), Element::NumDim());
  for (auto &source : Element::Sources()) {
    RealVec2 pnt (source->LocR(), source->LocS());
    s.col(0) += (source->fire(time) * Element::getDeltaFunctionCoefficients(pnt));
  }
  s.col(0) = Element::applyTestAndIntegrate(s.col(0));
  return s;
}

#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>
template class Elastic2D<TensorQuad<QuadP1>>;

