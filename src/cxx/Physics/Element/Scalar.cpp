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
RealMat Scalar<Element>::computeStress(const Ref<const RealMat> &strain) {

  // Interpolate the (square) of the velocity at each integration point.
  mVpSquared = Element::ParAtIntPts("VP").array().pow(2);

  // Calculate sigma_ux and sigma_uy.
  mStress.col(0) = mVpSquared.array().cwiseProduct(strain.col(0).array());
  mStress.col(1) = mVpSquared.array().cwiseProduct(strain.col(1).array());
  if (Element::NumDim() == 3) {
    mStress.col(2) = mVpSquared.array().cwiseProduct(strain.col(2).array());
  }
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
    RealVec pnt;
    if (Element::NumDim() == 2) { pnt.resize(2); pnt << source->LocR(), source->LocS(); }
    if (Element::NumDim() == 3) { pnt.resize(3); pnt << source->LocR(), source->LocS(), source->LocT(); }
    mSource += (source->fire(time) * Element::getDeltaFunctionCoefficients(pnt));
  }
  return Element::applyTestAndIntegrate(mSource);
}

#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/Hexahedra.h>
#include <Element/Simplex/Triangle.h>
#include <Element/Simplex/Tetrahedra.h>
#include <Element/HyperCube/HexP1.h>
#include <Element/HyperCube/QuadP1.h>
#include <Element/Simplex/TriP1.h>
#include <Element/Simplex/TetP1.h>

template class Scalar<TensorQuad<QuadP1>>;
template class Scalar<Hexahedra<HexP1>>;
template class Scalar<Triangle<TriP1>>;
template class Scalar<Tetrahedra<TetP1>>;
