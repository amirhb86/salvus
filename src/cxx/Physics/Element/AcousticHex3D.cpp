#include <Physics/AcousticHex3D.h>

// Dependencies.
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>
#include <Receiver/Receiver.h>
#include <Source/Source.h>
#include <Mesh/Mesh.h>

#include <iostream>

using namespace Eigen;

template <typename Element>
AcousticHex3D<Element>::AcousticHex3D(Options options): Element(options) {

  // Allocate all work arrays.
  mVpSquared.setZero(Element::NumIntPnt());
  mStiff.setZero(Element::NumIntPnt());
  mSource.setZero(Element::NumIntPnt());
  mStress.setZero(Element::NumIntPnt(), Element::NumDim());
  mStrain.setZero(Element::NumIntPnt(), Element::NumDim());

}

template <typename Element>
void AcousticHex3D<Element>::attachMaterialPropertiesNew(const ExodusModel *model) {
  Element::attachMaterialProperties(model, "VP");
}

template <typename Element>
std::vector<std::string> AcousticHex3D<Element>::PullElementalFields() const { return { "u" }; }

template <typename Element>
std::vector<std::string> AcousticHex3D<Element>::PushElementalFields() const { return { "a" }; }

template <typename Element>
void AcousticHex3D<Element>::assembleElementMassMatrix(Mesh *mesh) {

  // In this acoustic formulation we just multiply shape functions together.
  VectorXd mass_matrix = Element::applyTestAndIntegrate(VectorXd::Ones(Element::NumIntPnt()));
  
  // Sum up into global DOFs.  
  mesh->addFieldFromElement("m", Element::ElmNum(), Element::ClsMap(), mass_matrix);

}

MatrixXd AcousticHex3D<Element>::computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd> &u) {
  return Eigen::MatrixXd::Zero(Element::NumIntPnt(), 1);
}

template <typename Element>
double AcousticHex3D<Element>::CFL_estimate() {
  double vpMax = Element::ParAtIntPts("VP").maxCoeff();
  return Element::CFL_constant() * Element::estimatedElementRadius() / vpMax;
}

template <typename Element>
MatrixXd AcousticHex3D<Element>::computeStress(const Ref<const MatrixXd> &strain) {

  // Calculate sigma_ux and sigma_uy.
  mStress.col(0) = mVpSquared.array() * strain.col(0).array();
  mStress.col(1) = mVpSquared.array() * strain.col(1).array();
  mStress.col(2) = mVpSquared.array() * strain.col(2).array();
  return mStress;

}

template <typename Element>
void AcousticHex3D<Element>::prepareStiffness() {
  Element::precomputeConstants();
  mVpSquared = Element::ParAtIntPts("VP").array().pow(2);
}

template <typename Element>
MatrixXd AcousticHex3D<Element>::computeStiffnessTerm(
    const MatrixXd &u) {
  
  // Calculate gradient from displacement.
  // mStrain = Element::computeGradient(u.col(0));
  
  // Stress from strain.
  // mStress = computeStress(mStrain);
  
  // Complete application of K->u.
  // mStiff = Element::applyGradTestAndIntegrate(mStress);
  
  mStiff = Element::computeStiffnessFull(u.col(0),mVpSquared);
  
  return mStiff;

}

template <typename Element>
MatrixXd AcousticHex3D<Element>::computeSourceTerm(const double time) {
  mSource.setZero();  
  for (auto source : Element::Sources()) {
    mSource += (source->fire(time) * Element::getDeltaFunctionCoefficients(source->ReferenceLocationR(),
                                                                           source->ReferenceLocationS(),
                                                                           source->ReferenceLocationT()));
  }
  return mSource;
}


template <typename Element>
void AcousticHex3D<Element>::setupEigenfunctionTest(Mesh *mesh, Options options) {

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
double AcousticHex3D<Element>::checkEigenfunctionTest(Mesh *mesh, Options options,
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
  double vp = Element::ParAtIntPts("VP").mean();
  double un_t = cos(M_PI/Lx*sqrt(3)*time*vp);
  VectorXd exact = un_t * un_xyz;
  double element_error = (exact - u).array().abs().maxCoeff();
  if (!Element::ElmNum()) {
  }

  return element_error;

}

#include <Element/HyperCube/Hexahedra.h>
#include <Element/HyperCube/HexP1.h>
template class AcousticHex3D<Hexahedra<HexP1>>;

