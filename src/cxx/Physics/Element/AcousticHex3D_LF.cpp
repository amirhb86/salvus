#include <Physics/AcousticHex3D_LF.h>

// Dependencies.
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>
#include <Receiver/Receiver.h>
#include <Source/Source.h>
#include <Mesh/Mesh.h>

#include <iostream>

using namespace Eigen;

template <typename Element>
AcousticHex3D_LF<Element>::AcousticHex3D_LF(std::unique_ptr<Options> const &options): Element(options) {

  // Allocate all work arrays.
  mVpSquared.setZero(Element::NumIntPnt());
  mStiff.setZero(Element::NumIntPnt());
  mSource.setZero(Element::NumIntPnt());
  mStress.setZero(Element::NumIntPnt(), Element::NumDim());
  mStrain.setZero(Element::NumIntPnt(), Element::NumDim());

}

template <typename Element>
void AcousticHex3D_LF<Element>::attachMaterialProperties(std::unique_ptr<ExodusModel> const &model) {
  Element::attachMaterialProperties(model, "VP");
}

template <typename Element>
std::vector<std::string> AcousticHex3D_LF<Element>::PullElementalFields() const { return { "u" }; }

template <typename Element>
std::vector<std::string> AcousticHex3D_LF<Element>::PushElementalFields() const { return { "a" }; }

template <typename Element>
void AcousticHex3D_LF<Element>::assembleElementMassMatrix(Mesh *mesh) {

  // In this acoustic formulation we just multiply shape functions together.
  VectorXd mass_matrix = Element::applyTestAndIntegrate(VectorXd::Ones(Element::NumIntPnt()));
  
  // Sum up into global DOFs.  
  mesh->addFieldFromElement("m", Element::ElmNum(), Element::ClsMap(), mass_matrix);

}

template <typename Element>
MatrixXd AcousticHex3D_LF<Element>::computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd> &u) {
  return Eigen::MatrixXd::Zero(Element::NumIntPnt(), 1);
}

template <typename Element>
double AcousticHex3D_LF<Element>::CFL_estimate() {
  double vpMax = Element::ParAtIntPts("VP").maxCoeff();
  return Element::CFL_constant() * Element::estimatedElementRadius() / vpMax;
}

template <typename Element>
MatrixXd AcousticHex3D_LF<Element>::computeStress(const Ref<const MatrixXd> &strain) {

  // Calculate sigma_ux and sigma_uy.
  mStress.col(0) = mVpSquared.array() * strain.col(0).array();
  mStress.col(1) = mVpSquared.array() * strain.col(1).array();
  mStress.col(2) = mVpSquared.array() * strain.col(2).array();
  return mStress;

}

template <typename Element>
void AcousticHex3D_LF<Element>::prepareStiffness() {
  Element::precomputeConstants();
  mVpSquared = Element::ParAtIntPts("VP").array().pow(2);
}

template <typename Element>
MatrixXd AcousticHex3D_LF<Element>::computeStiffnessTerm(
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
MatrixXd AcousticHex3D_LF<Element>::computeSourceTerm(const double time) {
  mSource.setZero();  
  for (auto source : Element::Sources()) {
    mSource += (source->fire(time) * Element::getDeltaFunctionCoefficients(source->ReferenceLocationR(),
                                                                           source->ReferenceLocationS(),
                                                                           source->ReferenceLocationT()));
  }
  return mSource;
}

VectorXd exactSolution(VectorXd& xn, VectorXd& yn, VectorXd& zn, double time,
                       double L,
                       double x0, double y0, double z0,
                       double vp) {

  double Lx, Ly, Lz;
  Lx = Ly = Lz = L;
  VectorXd un_xyz =
    (M_PI/Lx*(xn.array()-(x0+L/2))).sin() *
    (M_PI/Ly*(yn.array()-(y0+L/2))).sin() *
    (M_PI/Lz*(zn.array()-(z0+L/2))).sin();

  double un_t = cos(M_PI/Lx*sqrt(3)*time*vp);

  return  un_t * un_xyz;
}

template <typename Element>
void AcousticHex3D_LF<Element>::setupEigenfunctionTest(Mesh *mesh, std::unique_ptr<Options> const &options) {

  double L, Lx, Ly, Lz;
  double x0 = options->IC_Center_x();
  double y0 = options->IC_Center_y();
  double z0 = options->IC_Center_z();
  
  L = Lx = Ly = Lz = options->IC_SquareSide_L();
  VectorXd pts_x, pts_y, pts_z;
  std::tie(pts_x, pts_y, pts_z) = Element::buildNodalPoints();
  // VectorXd un_xyz =
  //   (M_PI/Lx*(pts_x.array()-(x0+L/2))).sin() *
  //   (M_PI/Ly*(pts_y.array()-(y0+L/2))).sin() * 
  //   (M_PI/Lz*(pts_z.array()-(z0+L/2))).sin();
  auto dt = options->TimeStep();
  // auto time = -dt;
  double vp = Element::ParAtIntPts("VP").mean();
  // VectorXd unm1 = un_xyz*cos(-M_PI/Lx*sqrt(3)*dt*vp);
  VectorXd un = exactSolution(pts_x,pts_y,pts_z,0.0, Lx,x0,y0,z0, vp);
  VectorXd unm1 = exactSolution(pts_x,pts_y,pts_z,-dt, Lx,x0,y0,z0, vp);
  
  mesh->setFieldFromElement("u", Element::ElmNum(), Element::ClsMap(), un);
  mesh->setFieldFromElement("unm1", Element::ElmNum(), Element::ClsMap(), unm1);  
}

template <typename Element>
double AcousticHex3D_LF<Element>::checkEigenfunctionTest(Mesh *mesh, std::unique_ptr<Options> const &options,
                                                  const Ref<const MatrixXd>& u, double time) {

  double L, Lx, Ly, Lz;
  double x0 = options->IC_Center_x();
  double y0 = options->IC_Center_y();
  double z0 = options->IC_Center_z();
  L = Lx = Ly = Lz = options->IC_SquareSide_L();
  VectorXd pts_x, pts_y, pts_z;
  std::tie(pts_x,pts_y,pts_z) = Element::buildNodalPoints();
  // VectorXd un_xyz =
  //   (M_PI/Lx*(pts_x.array()-(x0+L/2))).sin() *
  //   (M_PI/Ly*(pts_y.array()-(y0+L/2))).sin() *
  //   (M_PI/Lz*(pts_z.array()-(z0+L/2))).sin();
  double vp = Element::ParAtIntPts("VP").mean();
  // double un_t = cos(M_PI/Lx*sqrt(3)*time*vp);
  // VectorXd exact = un_t * un_xyz;

  VectorXd exact = exactSolution(pts_x,pts_y,pts_z,time, Lx, x0,y0,z0, vp);
  
  mesh->setFieldFromElement("u_exact", Element::ElmNum(), Element::ClsMap(), exact);
  double element_error = (exact - u).array().abs().maxCoeff();
  if (!Element::ElmNum()) {
  }

  return element_error;

}

#include <Element/HyperCube/Hexahedra.h>
#include <Element/HyperCube/HexP1.h>
template class AcousticHex3D_LF<Hexahedra<HexP1>>;

