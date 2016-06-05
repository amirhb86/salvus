#include <Mesh/Mesh.h>
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>
#include <Physics/Elastic2D.h>
#include <Source/Source.h>

using namespace Eigen;

template <typename Element>
Elastic2D<Element>::Elastic2D(Options options): Element(options) {

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
void Elastic2D<Element>::attachMaterialPropertiesNew(const ExodusModel *model) {
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
void Elastic2D<Element>::assembleElementMassMatrix(Mesh *mesh) {

  VectorXd mass_matrix = Element::applyTestAndIntegrate(Element::ParAtIntPts("RHO"));
  mesh->addFieldFromElement("m", Element::ElmNum(), Element::ClsMap(), mass_matrix);

}

template <typename Element>
MatrixXd Elastic2D<Element>::computeStiffnessTerm(const Eigen::MatrixXd &u) {

  // strain ux_x, ux_y, uy_x, uy_y.
  mStrain.leftCols(2) = Element::computeGradient(u.col(0));
  mStrain.rightCols(2) = Element::computeGradient(u.col(1));

  /* TODO: PROPER CONVENTION!!!! */
  mc11 = Element::ParAtIntPts("C11");
  mc12 = Element::ParAtIntPts("C13");
//  mc13 = Element::ParAtIntPts("C12");
  mc22 = Element::ParAtIntPts("C33");
//  mc23 = Element::ParAtIntPts("C23");
  mc33 = Element::ParAtIntPts("C55");

  Eigen::VectorXd test = mStrain.col(1) + mStrain.col(2);

  mStress.col(0) =
      mc11.array().cwiseProduct(mStrain.col(0).array()) +
      mc12.array().cwiseProduct(mStrain.col(3).array()) +
      mc13.array().cwiseProduct(test.array());

  mStress.col(1) =
      mc12.array().cwiseProduct(mStrain.col(0).array()) +
      mc22.array().cwiseProduct(mStrain.col(3).array()) +
      mc23.array().cwiseProduct(test.array());

  mStress.col(2) =
      mc13.array().cwiseProduct(mStrain.col(0).array()) +
      mc23.array().cwiseProduct(mStrain.col(3).array()) +
      mc33.array().cwiseProduct(test.array());


  Eigen::MatrixXd temp(mStress.rows(), 2);
  temp.col(0) = mStress.col(0);
  temp.col(1) = mStress.col(2);

  mStiff.col(0) = Element::applyGradTestAndIntegrate(temp);

  temp.col(0) = mStress.col(2);
  temp.col(1) = mStress.col(1);
  mStiff.col(1) = Element::applyGradTestAndIntegrate(temp);


  return mStiff;

}

template <typename Element>
MatrixXd Elastic2D<Element>::computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd> &u) {
  return Eigen::MatrixXd::Zero(Element::NumIntPnt(), Element::NumDim());
}

template <typename Element>
MatrixXd Elastic2D<Element>::computeSourceTerm(const double time) {
  MatrixXd s = MatrixXd::Zero(Element::NumIntPnt(), Element::NumDim());
  for (auto &source : Element::Sources()) {
    s.col(0) += (source->fire(time) * Element::getDeltaFunctionCoefficients(
        source->ReferenceLocationR(), source->ReferenceLocationS()));
  }
  return s;
}

template <typename Element>
void Elastic2D<Element>::setupEigenfunctionTest(Mesh *mesh, Options options) {
  
  
  double L, Lx, Ly, Lz;
  double x0 = options.IC_Center_x();
  double y0 = options.IC_Center_y();
  double z0 = options.IC_Center_z();
  
  L = Lx = Ly = Lz = options.IC_SquareSide_L();
  VectorXd pts_x, pts_y;
  std::tie(pts_x, pts_y) = Element::buildNodalPoints();

  // # P wave
  // #un[:,1] = sin(pi / 2 * x_n)
  // #unm1[:,1] = cos(pi / 2 * sqrt(3) * dt) * sin(pi / 2 * x_n)
  // Eigen::VectorXd un_z(pts_x.size()); un_z.setZero();
  // Eigen::VectorXd un_x = (PI/Lx*(pts_x.array()-(x0+L/2))).sin();
    
  // # S wave in z direction (from Julia)
  // un[:,2] = sin(pi / 2 * x_n)
  // unm1[:,2] = cos(pi / 2 * dt) * sin(pi / 2 * x_n)    
  Eigen::VectorXd un_x(pts_x.size()); un_x.setZero();
  Eigen::VectorXd un_z = (M_PI/Lx*(pts_x.array()-(x0+L/2))).sin();

  mesh->setFieldFromElement("ux", Element::ElmNum(), Element::ClsMap(), un_x);    
  mesh->setFieldFromElement("uy", Element::ElmNum(), Element::ClsMap(), un_z);
  
};

template <typename Element>
double Elastic2D<Element>::checkEigenfunctionTest(Mesh *mesh, Options options,
                                                    const Ref<const MatrixXd>& u, double time) {


  double L, Lx, Ly, Lz;
  double x0 = options.IC_Center_x();
  double y0 = options.IC_Center_y();
  double z0 = options.IC_Center_z();
  L = Lx = Ly = Lz = options.IC_SquareSide_L();
  VectorXd pts_x, pts_y, pts_z;
  std::tie(pts_x,pts_y) = Element::buildNodalPoints();
  Eigen::MatrixXd un_xz(pts_x.size(),2);
  un_xz.col(0).setZero();
  un_xz.col(1) = (M_PI/Lx*(pts_x.array()-(x0+L/2))).sin();    
      
  double c11 = Element::ParAtIntPts("C11").mean();
  double c13 = Element::ParAtIntPts("C13").mean();
  double c15 = Element::ParAtIntPts("C15").mean();
  double c33 = Element::ParAtIntPts("C33").mean();
  double c35 = Element::ParAtIntPts("C35").mean();
  double c55 = Element::ParAtIntPts("C55").mean();
  double rho = Element::ParAtIntPts("Rho").mean();

  double VP = sqrt(c11/rho);
  double VS = sqrt(c33/rho);

  // s-wave only
  double un_t = cos(M_PI/Lx*time*VS);
  return (u - un_t*un_xz).array().abs().maxCoeff();
  
};

#include <Element/HyperCube/Quad.h>
#include <Element/HyperCube/QuadP1.h>
template class Elastic2D<Quad<QuadP1>>;

