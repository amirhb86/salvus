#include <Mesh/Mesh.h>
#include <Source/Source.h>
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>
#include <Physics/AcousticTet.h>

using namespace Eigen;

template <typename Element>
AcousticTet<Element>::AcousticTet(Options options): Element(options) {

  // Allocate all work arrays.
  mVpSquared.setZero(Element::NumIntPnt());
  mStiff.setZero(Element::NumIntPnt());
  mSource.setZero(Element::NumIntPnt());
  mStress.setZero(Element::NumIntPnt(), Element::NumDim());
  mStrain.setZero(Element::NumIntPnt(), Element::NumDim());

}

template <typename Element>
void AcousticTet<Element>::attachMaterialPropertiesNew(const ExodusModel *model) {
  Element::attachMaterialProperties(model, "VP");
}

template <typename Element>
std::vector<std::string> AcousticTet<Element>::PullElementalFields() const { return { "u" }; }

template <typename Element>
std::vector<std::string> AcousticTet<Element>::PushElementalFields() const { return { "a" }; }

template <typename Element>
void AcousticTet<Element>::assembleElementMassMatrix(Mesh *mesh) {

  // In this acoustic formulation we just multiply shape functions together.
  VectorXd mass_matrix = Element::applyTestAndIntegrate(VectorXd::Ones(Element::NumIntPnt()));
  
  // Sum up into global DOFs.  
  mesh->addFieldFromElement("m", Element::ElmNum(), Element::ClsMap(), mass_matrix);

}


template <typename Element>
MatrixXd AcousticTet<Element>::computeStress(const Ref<const MatrixXd> &strain) {

  // Interpolate the velocity at each integration point.
  auto vp = Element::ParAtIntPts("VP");

  // Calculate sigma_ux and sigma_uy. (note vp^2)
  mStress.col(0) = vp.array().pow(2).cwiseProduct(strain.col(0).array());
  mStress.col(1) = vp.array().pow(2).cwiseProduct(strain.col(1).array());
  mStress.col(2) = vp.array().pow(2).cwiseProduct(strain.col(2).array());
  
  // 5% faster
  // for(int i=0;i<strain.rows();i++) {
  //   double vpi = vp(i);
  //   mStress(i,0) *= vpi*vpi;
  //   mStress(i,1) *= vpi*vpi;
  //   mStress(i,2) *= vpi*vpi;
  // }

  return mStress;

}

template <typename Element>
MatrixXd AcousticTet<Element>::computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd> &u) {
  return Eigen::MatrixXd::Zero(Element::NumIntPnt(), 1);
}


template <typename Element>
void AcousticTet<Element>::prepareStiffness() {
  auto velocity = Element::ParAtIntPts("VP");
  mElementStiffnessMatrix = Element::buildStiffnessMatrix(velocity);
}

template <typename Element>
MatrixXd AcousticTet<Element>::computeStiffnessTerm(
    const MatrixXd &u) {

  // Ku version
  mStiff.noalias() = mElementStiffnessMatrix*u.col(0);
  
  // Calculate gradient from displacement.
  mStrain = Element::computeGradient(u.col(0));
  
  // Seperated version
  // Stress from strain.
  // mStress = computeStress(mStrain);
  
  // Complete application of K->u.
  // mStiff = Element::applyGradTestAndIntegrate(mStress);
  
  // if(Element::ElmNum() == 0) {
  //   std::cout << "mStiffv1=" << mStiff.transpose() << "\n";
  //   std::cout << "mStiffv2=" << mStiffv2.transpose() << "\n";
  // }
  //   VectorXd delta; delta.setZero(mStiffv2.size());
  //   delta(0) = 1.0;
  //   MatrixXd mStrain_delta = Element::computeGradient(delta);
  //   MatrixXd mStress_delta = computeStress(mStrain_delta);
  //   VectorXd mStiff_delta = Element::applyGradTestAndIntegrate(mStress_delta);
  //   std::cout << "K[:,0] = " << mElementStiffnessMatrix.col(0).transpose() << "\n";
  //   std::cout << "mStiff_delta=" << mStiff_delta.transpose() << "\n";
  // }
  
  return mStiff;
}

template <typename Element>
MatrixXd AcousticTet<Element>::computeSourceTerm(const double time) {
  mSource.setZero();  
  for (auto source : Element::Sources()) {
    mSource += (source->fire(time) * Element::getDeltaFunctionCoefficients(source->ReferenceLocationR(),
                                                                           source->ReferenceLocationS(),
                                                                           source->ReferenceLocationT()));
  }
  return mSource;
}


template <typename Element>
void AcousticTet<Element>::setupEigenfunctionTest(Mesh *mesh, Options options) {

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
double AcousticTet<Element>::checkEigenfunctionTest(Mesh *mesh, Options options,
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

#include <Element/Simplex/Tetrahedra.h>
#include <Element/Simplex/TetP1.h>
template class AcousticTet<Tetrahedra<TetP1>>;

