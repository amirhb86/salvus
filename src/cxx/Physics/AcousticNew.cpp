#include "AcousticNew.h"

#include <Mesh/Mesh.h>
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>

using namespace Eigen;

template <typename Element>
AcousticNew<Element>::AcousticNew(Options options): Element(options) {

  // Allocate all work arrays.
  mVpSquared.setZero(Element::NumIntPnt());
  mStiff.setZero(Element::NumIntPnt());
  mStress.setZero(Element::NumIntPnt(), Element::NumDim());
  mStrain.setZero(Element::NumIntPnt(), Element::NumDim());

}

template <typename Element>
void AcousticNew<Element>::attachMaterialPropertiesNew(const ExodusModel *model) {
  Element::attachMaterialProperties(model, "VP");
}

template <typename Element>
std::vector<std::string> AcousticNew<Element>::PullElementalFields() const { return { "u" }; }

template <typename Element>
std::vector<std::string> AcousticNew<Element>::PushElementalFields() const { return { "a" }; }

template <typename Element>
void AcousticNew<Element>::assembleElementMassMatrix(Mesh *mesh) {

  // In this acoustic formulation we just multiply shape functions together.
  VectorXd mass_matrix = Element::applyTestAndIntegrate(VectorXd::Ones(Element::NumIntPnt()));

  // Sum up into global DOFs.
  mesh->addFieldFromElement("m", Element::ElmNum(), Element::ClsMap(), mass_matrix);

}

template <typename Element>
MatrixXd AcousticNew<Element>::computeStiffnessTerm(
    const MatrixXd &u) {

  // Calculate gradient from displacement.
  mStrain = Element::computeGradient(u);

  // Interpolate the (square) of velocity to each integration point.
  mVpSquared = Element::ParAtIntPts("VP").array().pow(2);

  // Calculate sigma_ux and sigma_uy.
  mStress.col(0) = mVpSquared.array().cwiseProduct(mStrain.col(0).array());
  mStress.col(1) = mVpSquared.array().cwiseProduct(mStrain.col(1).array());

  // Complete application of K->u.
  mStiff = Element::applyGradTestAndIntegrate(mStress);

  return mStiff;

}

template <typename Element>
MatrixXd AcousticNew<Element>::computeSourceTerm(const double time) { return Eigen::MatrixXd(1,1); }

template <typename Element>
void AcousticNew<Element>::setupEigenfunctionTest(Mesh *mesh, Options options) {

  double L, Lx, Ly;
  double x0 = options.IC_Center_x();
  double y0 = options.IC_Center_z();
  L = Lx = Ly = options.IC_SquareSide_L();
  VectorXd pts_x, pts_y;
  std::tie(pts_x, pts_y) = Element::buildNodalPoints();
  VectorXd un = (M_PI/Lx*(pts_x.array()-(x0+L/2))).sin() * (M_PI/Ly*(pts_y.array()-(y0+L/2))).sin();
  VectorXd vn = VectorXd::Zero(pts_x.size());
  VectorXd an = VectorXd::Zero(pts_x.size());
  mesh->setFieldFromElement("u", Element::ElmNum(), Element::ClsMap(), un);
  mesh->setFieldFromElement("v", Element::ElmNum(), Element::ClsMap(), vn);
  mesh->setFieldFromElement("a_", Element::ElmNum(), Element::ClsMap(), an);

}

template <typename Element>
double AcousticNew<Element>::checkEigenfunctionTest(Mesh *mesh, Options options,
                                                  const Ref<const MatrixXd>& u, double time) {

  double L, Lx, Ly;
  double x0 = options.IC_Center_x();
  double y0 = options.IC_Center_z();
  L = Lx = Ly = options.IC_SquareSide_L();
  VectorXd pts_x, pts_y;
  std::tie(pts_x,pts_y) = Element::buildNodalPoints();
  VectorXd un_xy = (M_PI/Lx*(pts_x.array()-(x0+L/2))).sin()*(M_PI/Ly*(pts_y.array()-(y0+L/2))).sin();
  double vp = Element::ParAtIntPts("VP").mean();
  double un_t = cos(M_PI/Lx*sqrt(2)*time*vp);
  VectorXd exact = un_t * un_xy;
  double element_error = (exact - u).array().abs().maxCoeff();
  if (!Element::ElmNum()) {
  }

  return element_error;

}

#include <QuadNew.h>
#include <QuadP1.h>
template class AcousticNew<QuadNew<QuadP1>>;

