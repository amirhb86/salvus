//
// Created by Michael Afanasiev on 2016-05-03.
//

#include <Element/HyperCube/Quad/QuadP1.h>
#include <Element/HyperCube/QuadNew.h>
#include "AcousticNew.h"

using namespace Eigen;

template <typename Derived>
AcousticNew<Derived>::AcousticNew(Options options): Derived(options) {

  mStiff.setZero(Derived::NumIntPnt());
  mVelGrd.setZero(Derived::NumIntPnt());
  mDetJac.setZero(Derived::NumIntPnt());
  mStress.setZero(Derived::NumIntPnt(), Derived::NumDimNew());
  mStrain.setZero(Derived::NumIntPnt(), Derived::NumDimNew());

}

template <typename Derived>
void AcousticNew<Derived>::attachMaterialPropertiesNew(ExodusModel *model) {
  Derived::attachMaterialPropertiesNew(model, "VP");
}

template <typename Derived>
std::vector<std::string> AcousticNew<Derived>::PullElementalFields() const { return { "u" }; }

template <typename Derived>
std::vector<std::string> AcousticNew<Derived>::PushElementalFields() const { return { "a" }; }

template <typename Derived>
void AcousticNew<Derived>::assembleElementMassMatrix(Mesh *mesh) {

  int i = 0;
  double detJac;
  Matrix2d invJac;
  VectorXd mass_matrix(Derived::NumIntPnt());
  for (int s_ind = 0; s_ind < Derived::NumIntPtsS(); s_ind++) {
    for (int r_ind = 0; r_ind < Derived::NumIntPtsR(); r_ind++) {

      // (r,s) coordinates for this point
      double r = Derived::IntCrdR(r_ind);
      double s = Derived::IntCrdS(s_ind);

      // mass value at point.
      std::tie(invJac, detJac) = Derived::inverseJacobianAtPoint(r, s, Derived::VtxCrd());
      mass_matrix(i) = detJac * Derived::IntWgtR(r_ind) * Derived::IntWgtS(s_ind);
      i++;

    }
  }

  // Sum up into global DOFs.
  mesh->addFieldFromElement("m", Derived::ElmNum(), Derived::ClosureMap(), mass_matrix);
}

template <typename Derived>
MatrixXd AcousticNew<Derived>::computeStiffnessTerm(
    const MatrixXd &displacement) {

  int itr = 0;
  int s_pts = Derived::NumIntPtsS();
  int r_pts = Derived::NumIntPtsR();
  Matrix2d invJac;
  Vector2d refStrain;

  for (int s_ind = 0; s_ind < s_pts; s_ind++) {
    for (int r_ind = 0; r_ind < r_pts; r_ind++) {

      // (r,s) coordinates for this point.
      double r = Derived::IntCrdR(r_ind);
      double s = Derived::IntCrdS(s_ind);

      // strain.
      std::tie(invJac, mDetJac(itr)) = Derived::inverseJacobianAtPoint(r, s, Derived::VtxCrd());
      mStrain.row(itr) = invJac * (refStrain <<
        Derived::GrdRow(r_ind).dot(Derived::rVectorStride(displacement, s_ind, s_pts, r_pts)),
        Derived::GrdRow(s_ind).dot(Derived::sVectorStride(displacement, r_ind, s_pts, r_pts))).finished();

      // stress.
      double vp = Derived::ParAtPnt(r, s, "VP");
      mStress.row(itr) = vp * vp * mStrain.row(itr);

      // advance.
      itr++;

    }
  }

  itr = 0;
  for (int s_ind = 0; s_ind < s_pts; s_ind++) {
    for (int r_ind = 0; r_ind < r_pts; r_ind++) {

      // (r,s) coordiantes for this point.
      double r = Derived::IntCrdR(r_ind);
      double s = Derived::IntCrdS(s_ind);

      // get inverse jacobian again.
      double _;
      std::tie(invJac, _) = Derived::inverseJacobianAtPoint(r, s, Derived::VtxCrd());

      double dphi_r_dux = Derived::IntWgtS(s_ind) *
          Derived::VecIntWgtR().dot(((Derived::rVectorStride(mDetJac, s_ind, s_pts, r_pts)).array() *
          Derived::rVectorStride(mStress.col(0), s_ind, s_pts, r_pts).array() *
          Derived::GrdCol(r_ind).array()).matrix());

      double dphi_s_dux = Derived::IntWgtR(r_ind) *
          Derived::VecIntWgtS().dot(((Derived::sVectorStride(mDetJac, r_ind, s_pts, r_pts)).array() *
          Derived::sVectorStride(mStress.col(0), r_ind, s_pts, r_pts).array() *
          Derived::GrdCol(s_ind).array()).matrix());

      double dphi_r_duy = Derived::IntWgtS(s_ind) *
          Derived::VecIntWgtR().dot(((Derived::rVectorStride(mDetJac, s_ind, s_pts, r_pts)).array() *
          Derived::rVectorStride(mStress.col(1), s_ind, s_pts, r_pts).array() *
          Derived::GrdCol(r_ind).array()).matrix());

      double dphi_s_duy = Derived::IntWgtR(r_ind) *
          Derived::VecIntWgtS().dot(((Derived::sVectorStride(mDetJac, r_ind, s_pts, r_pts)).array() *
          Derived::sVectorStride(mStress.col(1), r_ind, s_pts, r_pts).array() *
          Derived::GrdCol(s_ind).array()).matrix());

      Vector2d dphi_epseta_dux, dphi_epseta_duy;
      dphi_epseta_dux << dphi_r_dux, dphi_s_dux;
      dphi_epseta_duy << dphi_s_dux, dphi_s_duy;

      mStiff(itr) =
          invJac.row(0).dot(dphi_epseta_dux) +
          invJac.row(1).dot(dphi_epseta_duy);

      itr++;

    }
  }

  return mStiff;

}

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
  mesh->setFieldFromElement("u", Element::ElmNum(), Element::ClosureMap(), un);
  mesh->setFieldFromElement("v", Element::ElmNum(), Element::ClosureMap(), vn);
  mesh->setFieldFromElement("a_", Element::ElmNum(), Element::ClosureMap(), an);
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
  double vp = Element::ParAtPnt(0, 0, "VP");
  double un_t = cos(M_PI/Lx*sqrt(2)*time*vp);
  VectorXd exact = un_t * un_xy;
  double element_error = (exact - u).array().abs().maxCoeff();
  if (!Element::ElmNum()) {
//    std::cout << "EXACT:\n" << exact << "\nU:\n" << u << std::endl;
  }

  return element_error;
}

template class AcousticNew<QuadNew<QuadP1>>;

