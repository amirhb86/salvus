#include <Mesh/Mesh.h>
#include <Source/Source.h>
#include <Physics/Scalar.h>
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>

using namespace Eigen;

template <typename Element>
Scalar<Element>::Scalar(std::unique_ptr<Options> const &options): Element(options) {

  // Allocate all work arrays.
  mVpSquared.setZero(Element::NumIntPnt());
  mStiff.setZero(Element::NumIntPnt());
  mSource.setZero(Element::NumIntPnt());
  mStress.setZero(Element::NumIntPnt(), Element::NumDim());
  mStrain.setZero(Element::NumIntPnt(), Element::NumDim());

}

template <typename Element>
void Scalar<Element>::attachMaterialProperties(std::unique_ptr<ExodusModel> const &model) {
  Element::attachMaterialProperties(model, "VP");
}

template <typename Element>
std::vector<std::string> Scalar<Element>::PullElementalFields() const { return { "u" }; }

template <typename Element>
std::vector<std::string> Scalar<Element>::PushElementalFields() const { return { "a" }; }

template <typename Element>
RealMat Scalar<Element>::assembleElementMassMatrix() {

  // In this acoustic formulation we just multiply shape functions together.
  return Element::applyTestAndIntegrate(VectorXd::Ones(Element::NumIntPnt()));

}

template <typename Element>
double Scalar<Element>::CFL_estimate() {
//  double vpMax = Element::ParAtIntPts("VP").maxCoeff();
//  return Element::CFL_constant() * Element::estimatedElementRadius() / vpMax;
  return 1.0;
}


template <typename Element>
RealMat Scalar<Element>::computeStress(const Ref<const RealMat> &strain) {

  // Interpolate the (square) of the velocity at each integration point.
  mVpSquared = Element::ParAtIntPts("VP").array().pow(2);

  // Calculate sigma_ux and sigma_uy.
  mStress.col(0) = mVpSquared.array().cwiseProduct(strain.col(0).array());
  mStress.col(1) = mVpSquared.array().cwiseProduct(strain.col(1).array());
  return mStress;

}

template <typename Element>
RealMat Scalar<Element>::computeStiffnessTerm(const Ref<const RealMat>& u) {

  // Calculate gradient from displacement.
  mStrain = Element::computeGradient(u.col(0));

  // Stress from strain.
  mStress = computeStress(mStrain);

  // Complete application of K->u.
  mStiff = Element::applyGradTestAndIntegrate(mStress);

  return mStiff;

}

template <typename Element>
RealMat Scalar<Element>::computeSurfaceIntegral(const Ref<const RealMat> &u) {
  return RealMat::Zero(Element::NumIntPnt(), 1);
}

template <typename Element>
MatrixXd Scalar<Element>::computeSourceTerm(const double time) {
  mSource.setZero();
  for (auto &source : Element::Sources()) {
    mSource += (source->fire(time) * Element::getDeltaFunctionCoefficients(
        source->LocR(), source->LocS()));
  }
  return mSource;
}


//template <typename Element>
//void Scalar<Element>::setupEigenfunctionTest(std::unique_ptr<Mesh> const &mesh, std::unique_ptr<Options> const &options) {
//
//  double L, Lx, Ly;
//  double x0 = options->IC_Center_x();
//  double y0 = options->IC_Center_z();
//  L = Lx = Ly = options->IC_SquareSide_L();
//  VectorXd pts_x, pts_y;
//  std::tie(pts_x, pts_y) = Element::buildNodalPoints();
//  VectorXd un = (M_PI/Lx*(pts_x.array()-(x0+L/2))).sin() * (M_PI/Ly*(pts_y.array()-(y0+L/2))).sin();
//  VectorXd vn = VectorXd::Zero(pts_x.size());
//  VectorXd an = VectorXd::Zero(pts_x.size());
//  mesh->setFieldFromElement("u", Element::ElmNum(), Element::ClsMap(), un);
//  mesh->setFieldFromElement("v", Element::ElmNum(), Element::ClsMap(), vn);
//  mesh->setFieldFromElement("a_", Element::ElmNum(), Element::ClsMap(), an);
//
//}
//
//template <typename Element>
//double Scalar<Element>::checkEigenfunctionTest(std::unique_ptr<Mesh> const &mesh, std::unique_ptr<Options> const &options,
//                                                  const Ref<const MatrixXd>& u, double time) {
//
//  double L, Lx, Ly;
//  double x0 = options->IC_Center_x();
//  double y0 = options->IC_Center_z();
//  L = Lx = Ly = options->IC_SquareSide_L();
//  VectorXd pts_x, pts_y;
//  std::tie(pts_x,pts_y) = Element::buildNodalPoints();
//  VectorXd un_xy = (M_PI/Lx*(pts_x.array()-(x0+L/2))).sin()*(M_PI/Ly*(pts_y.array()-(y0+L/2))).sin();
//  double vp = Element::ParAtIntPts("VP").mean();
//  double un_t = cos(M_PI/Lx*sqrt(2)*time*vp);
//  VectorXd exact = un_t * un_xy;
//  double element_error = (exact - u).array().abs().maxCoeff();
//  if (!Element::ElmNum()) {
//  }
//
//  return element_error;
//
//}

#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>
template class Scalar<TensorQuad<QuadP1>>;
