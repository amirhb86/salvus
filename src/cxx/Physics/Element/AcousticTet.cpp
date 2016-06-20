#include <Mesh/Mesh.h>
#include <Source/Source.h>
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>
#include <Physics/AcousticTet.h>

using namespace Eigen;

template <typename Element>
AcousticTet<Element>::AcousticTet(std::unique_ptr<Options> const &options): Element(options) {

  // Allocate all work arrays.
  mVpSquared.setZero(Element::NumIntPnt());
  mStiff.setZero(Element::NumIntPnt());
  mSource.setZero(Element::NumIntPnt());
  mStress.setZero(Element::NumIntPnt(), Element::NumDim());
  mStrain.setZero(Element::NumIntPnt(), Element::NumDim());

}

template <typename Element>
void AcousticTet<Element>::attachMaterialProperties(std::unique_ptr<ExodusModel> const &model) {
  Element::attachMaterialProperties(model, "VP");
}

template <typename Element>
std::vector<std::string> AcousticTet<Element>::PullElementalFields() const { return { "u" }; }

template <typename Element>
std::vector<std::string> AcousticTet<Element>::PushElementalFields() const { return { "a" }; }

template <typename Element>
void AcousticTet<Element>::assembleElementMassMatrix(std::unique_ptr<Mesh> const &mesh) {

  // In this acoustic formulation we just multiply shape functions together.
  VectorXd mass_matrix = Element::applyTestAndIntegrate(VectorXd::Ones(Element::NumIntPnt()));
  
  // Sum up into global DOFs.  
  mesh->addFieldFromElement("m", Element::ElmNum(), Element::ClsMap(), mass_matrix);

}

template <typename Element>
double AcousticTet<Element>::CFL_estimate() {
  double vpMax = Element::ParAtIntPts("VP").maxCoeff();
  return Element::CFL_constant() * Element::estimatedElementRadius() / vpMax;
}


template <typename Element>
MatrixXd AcousticTet<Element>::computeStress(const Ref<const MatrixXd> &strain) {

  // Interpolate the velocity at each integration point.
  // auto vp = Element::ParAtIntPts("VP");
  
  // Calculate sigma_ux and sigma_uy. (note vp^2)
  mStress.col(0) = mVpSquared.cwiseProduct(strain.col(0));
  mStress.col(1) = mVpSquared.cwiseProduct(strain.col(1));
  mStress.col(2) = mVpSquared.cwiseProduct(strain.col(2));
  
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
  mVpSquared = velocity.array().pow(2);
  mElementStiffnessMatrix = Element::buildStiffnessMatrix(mVpSquared);
  
}

template <typename Element>
MatrixXd AcousticTet<Element>::computeStiffnessTerm(
    const MatrixXd &u) {

  // Ku version
  
  mStiff.noalias() = mElementStiffnessMatrix*u.col(0);
  
  // Calculate gradient from displacement (if needed)
  // mStrain = Element::computeGradient(u.col(0));
  
  // full stiffness "all-at-once"
  // mStiff = Element::computeStiffnessFull(u.col(0),mVpSquared);
  
  
  // // Seperated version
  // // Calculate gradient from displacement.
  // mStrain = Element::computeGradient(u.col(0));

  // // Stress from strain.
  // mStress = computeStress(mStrain);
  
  // Complete application of K->u.
  // mStiff = Element::applyGradTestAndIntegrate(mStress);

  return mStiff;
}

template <typename Element>
MatrixXd AcousticTet<Element>::computeSourceTerm(const double time) {
  mSource.setZero();  
  for (auto source : Element::Sources()) {
    mSource += (source->fire(time) * Element::getDeltaFunctionCoefficients(source->LocR(),
                                                                           source->LocS(),
                                                                           source->LocT()));
  }
  return mSource;
}


template <typename Element>
void AcousticTet<Element>::setupEigenfunctionTest(std::unique_ptr<Mesh> const &mesh, std::unique_ptr<Options> const &options) {

  double L, Lx, Ly, Lz;
  double x0 = options->IC_Center_x();
  double y0 = options->IC_Center_y();
  double z0 = options->IC_Center_z();
  
  L = Lx = Ly = Lz = options->IC_SquareSide_L();
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
double AcousticTet<Element>::checkEigenfunctionTest(std::unique_ptr<Mesh> const &mesh, std::unique_ptr<Options> const &options,
                                                  const Ref<const MatrixXd>& u, double time) {

  double L, Lx, Ly, Lz;
  double x0 = options->IC_Center_x();
  double y0 = options->IC_Center_y();
  double z0 = options->IC_Center_z();
  L = Lx = Ly = Lz = options->IC_SquareSide_L();
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

