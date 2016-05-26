#include "Acoustic3D_V.h"

// Dependencies.
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>
#include <Receiver/Receiver.h>
#include <Source/Source.h>
#include <Mesh/Mesh.h>

using namespace Eigen;

template <typename Element>
Acoustic3D_V<Element>::Acoustic3D_V(Options options): Element(options) {

  // Allocate all work arrays.
  mVpSquared.setZero(Element::NumIntPnt());
  mStiff.setZero(Element::NumIntPnt());
  mSource.setZero(Element::NumIntPnt());
  mStress.setZero(Element::NumIntPnt(), Element::NumDim());
  mStrain.setZero(Element::NumIntPnt(), Element::NumDim());

}

template <typename Element>
void Acoustic3D_V<Element>::attachMaterialPropertiesNew(const ExodusModel *model) {
  Element::attachMaterialProperties(model, "VPV");
}

template <typename Element>
std::vector<std::string> Acoustic3D_V<Element>::PullElementalFields() const { return { "u" }; }

template <typename Element>
std::vector<std::string> Acoustic3D_V<Element>::PushElementalFields() const { return { "a" }; }

template <typename Element>
void Acoustic3D_V<Element>::assembleElementMassMatrix(Mesh *mesh) {

  // In this acoustic formulation we just multiply shape functions together.
  VectorXd mass_matrix = Element::applyTestAndIntegrate(VectorXd::Ones(Element::NumIntPnt()));
  // std::cout << "clsmap:" << Element::ClsMap().transpose() << "\n";
  // Sum up into global DOFs.
  mesh->addFieldFromElement("m", Element::ElmNum(), Element::ClsMap(), mass_matrix);

}

template <typename Element>
MatrixXd Acoustic3D_V<Element>::computeStress(const Ref<const MatrixXd> &strain) {

  // Interpolate the (square) of the velocity at each integration point.
  mVpSquared = Element::ParAtIntPts("VPV").array().pow(2);

  // Calculate sigma_ux and sigma_uy.
  mStress.col(0) = mVpSquared.array().cwiseProduct(strain.col(0).array());
  mStress.col(1) = mVpSquared.array().cwiseProduct(strain.col(1).array());
  mStress.col(2) = mVpSquared.array().cwiseProduct(strain.col(2).array());
  return mStress;

}

template <typename Element>
MatrixXd Acoustic3D_V<Element>::computeStiffnessTerm(
    const MatrixXd &u) {

  // Calculate gradient from displacement.
  mStrain = Element::computeGradient(u);

  // Stress from strain.
  mStress = computeStress(mStrain);

  // Complete application of K->u.
  mStiff = Element::applyGradTestAndIntegrate(mStress);

  return mStiff;

}

template <typename Element>
MatrixXd Acoustic3D_V<Element>::computeSourceTerm(const double time) {
  mSource.setZero();  
  for (auto source : Element::Sources()) {
    mSource += (source->fire(time) * Element::getDeltaFunctionCoefficients(source->ReferenceLocationR(),
                                                                           source->ReferenceLocationS(),
                                                                           source->ReferenceLocationT()));
  }
  return mSource;
}


template <typename Element>
void Acoustic3D_V<Element>::setupEigenfunctionTest(Mesh *mesh, Options options) {

  double L, Lx, Ly, Lz;
  double x0 = options.IC_Center_x();
  double y0 = options.IC_Center_y();
  double z0 = options.IC_Center_z();
  
  L = Lx = Ly = Lz = options.IC_SquareSide_L();
  VectorXd pts_x, pts_y, pts_z;
  std::tie(pts_x, pts_y, pts_z) = Element::buildNodalPoints();
  VectorXd un =
    (M_PI/Lx*(pts_x.array()-(x0+L/2))).sin() *
    (M_PI/Ly*(pts_y.array()-(y0+L/2))).sin() * 
    (M_PI/Lz*(pts_z.array()-(z0+L/2))).sin();
  VectorXd vn = VectorXd::Zero(pts_x.size());
  VectorXd an = VectorXd::Zero(pts_x.size());
  mesh->setFieldFromElement("u", Element::ElmNum(), Element::ClsMap(), un);
  mesh->setFieldFromElement("v", Element::ElmNum(), Element::ClsMap(), vn);
  mesh->setFieldFromElement("a_", Element::ElmNum(), Element::ClsMap(), an);

}

template <typename Element>
double Acoustic3D_V<Element>::checkEigenfunctionTest(Mesh *mesh, Options options,
                                                  const Ref<const MatrixXd>& u, double time) {

  double L, Lx, Ly, Lz;
  double x0 = options.IC_Center_x();
  double y0 = options.IC_Center_y();
  double z0 = options.IC_Center_z();
  L = Lx = Ly = Lz = options.IC_SquareSide_L();
  VectorXd pts_x, pts_y, pts_z;
  std::tie(pts_x,pts_y,pts_z) = Element::buildNodalPoints();
  VectorXd un_xyz =
    (M_PI/Lx*(pts_x.array()-(x0+L/2))).sin() *
    (M_PI/Ly*(pts_y.array()-(y0+L/2))).sin() *
    (M_PI/Lz*(pts_z.array()-(z0+L/2))).sin();
  double vp = Element::ParAtIntPts("VPV").mean();
  double un_t = cos(M_PI/Lx*sqrt(3)*time*vp);
  VectorXd exact = un_t * un_xyz;
  double element_error = (exact - u).array().abs().maxCoeff();
  if (!Element::ElmNum()) {
  }

  return element_error;

}

#include <Hexahedra.h>
#include <HexP1.h>
template class Acoustic3D_V<Hexahedra<HexP1>>;
