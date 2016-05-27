#include "ElasticNew.h"

#include <Mesh/Mesh.h>
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>

using namespace Eigen;

template <typename Element>
ElasticNew<Element>::ElasticNew(Options options): Element(options) {

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
  mStrain.setZero(Element::NumIntPnt(), num_grad_cmps);

}

template <typename Element>
void ElasticNew<Element>::attachMaterialPropertiesNew(const ExodusModel *model) {
  Element::attachMaterialProperties(model, "RHO");
  Element::attachMaterialProperties(model, "C11");
  Element::attachMaterialProperties(model, "C12");
  Element::attachMaterialProperties(model, "C13");
  Element::attachMaterialProperties(model, "C22");
  Element::attachMaterialProperties(model, "C23");
  Element::attachMaterialProperties(model, "C33");
}

template <typename Element>
std::vector<std::string> ElasticNew<Element>::PullElementalFields() const { return { "ux", "uy"}; }

template <typename Element>
std::vector<std::string> ElasticNew<Element>::PushElementalFields() const { return { "ax", "ay"}; }

template <typename Element>
void ElasticNew<Element>::assembleElementMassMatrix(Mesh *mesh) {

  VectorXd mass_matrix = Element::applyTestAndIntegrate(Element::ParAtIntPts("RHO"));
  mesh->addFieldFromElement("m", Element::ElmNum(), Element::ClsMap(), mass_matrix);

}

template <typename Element>
MatrixXd ElasticNew<Element>::computeStiffnessTerm(const Eigen::MatrixXd &u) {

  // strain ux_x, uy_y, uy_x.
  mStrain.col(0) = Element::computeGradient(u.col(0)).col(0);
  mStrain.rightCols<2>() = Element::computeGradient(u.col(1));

  mc11 = Element::ParAtIntPts("C11");
  mc12 = Element::ParAtIntPts("C12");
  mc13 = Element::ParAtIntPts("C13");
  mc22 = Element::ParAtIntPts("C22");
  mc23 = Element::ParAtIntPts("C23");
  mc33 = Element::ParAtIntPts("C33");

  mStress.col(0) =
      mc11.array().cwiseProduct(mStrain.col(0).array()) +
      mc12.array().cwiseProduct(mStrain.col(1).array()) +
      mc13.array().cwiseProduct(2*mStrain.col(2).array());

  mStress.col(1) =
      mc12.array().cwiseProduct(mStrain.col(0).array()) +
      mc22.array().cwiseProduct(mStrain.col(1).array()) +
      mc23.array().cwiseProduct(2*mStrain.col(2).array());

  mStress.col(2) =
      mc13.array().cwiseProduct(mStrain.col(0).array()) +
      mc23.array().cwiseProduct(mStrain.col(1).array()) +
      mc33.array().cwiseProduct(2*mStrain.col(2).array());

  mStiff.col(0) = Element::applyGradTestAndIntegrate(
      Map<VectorXd,0,OuterStride<>>(mStress.data(), Element::NumIntPnt(),
                                    OuterStride<>(Element::NumIntPnt())));
  mStiff.col(1) = Element::applyGradTestAndIntegrate(
      mStress.rightCols<2>());

  return mStiff;

}

template <typename Element>
MatrixXd ElasticNew<Element>::computeSourceTerm(const double time) { return Eigen::MatrixXd(1,1); }

template <typename Element>
void ElasticNew<Element>::setupEigenfunctionTest(Mesh *mesh, Options options) {};

template <typename Element>
double ElasticNew<Element>::checkEigenfunctionTest(Mesh *mesh, Options options,
                                                    const Ref<const MatrixXd>& u, double time) { return 0.0; };

#include <QuadNew.h>
#include <QuadP1.h>
template class ElasticNew<QuadNew<QuadP1>>;

