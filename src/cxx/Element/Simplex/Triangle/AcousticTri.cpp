#include "AcousticTri.h"

AcousticTri::AcousticTri(Options options): Triangle(options) {
    
    // Allocate element vectors.
    mMssMat.setZero(mNumIntPnt);

    // Strain matrix.
    mElementStrain.setZero(2, mNumIntPnt);

}

void AcousticTri::buildStiffnessMatrix() {

  Eigen::Matrix2d invJ;
  double detJ;
  std::tie(invJ,detJ) = inverseJacobianAtPoint(0,0);
  //Jinv= rx, sx,
  //      rz, sz;
  auto drdx = invJ(0,0);
  auto dsdx = invJ(0,1);
  auto drdz = invJ(1,0);
  auto dsdz = invJ(1,1);

  // build material on all nodes
  Eigen::VectorXd velocity(mNumIntPnt);
  mElementStiffnessMatrix.resize(mNumIntPnt,mNumIntPnt);
    
  // interpolate velocity at all nodes
  for(int i=0;i<mNumIntPnt;i++) {
    auto r = mIntegrationCoordinates_r[i];
    auto s = mIntegrationCoordinates_s[i];
    velocity(i) = interpolateAtPoint(r,s).dot(mMaterialVelocityAtVertices);
  }
    
  // loop over matrix(i,j)
  for(int i=0;i<mNumIntPnt;i++) {
      
    Eigen::VectorXd dPhi_dr_i = mGradientPhi_dr.row(i);
    Eigen::VectorXd dPhi_ds_i = mGradientPhi_ds.row(i);
    auto dPhi_dx_i = dPhi_dr_i*drdx + dPhi_ds_i*dsdx;
    auto dPhi_dz_i = dPhi_dr_i*drdz + dPhi_ds_i*dsdz;
    for(int j=0;j<mNumIntPnt;j++) {
      Eigen::VectorXd dPhi_dr_j = mGradientPhi_dr.row(j);
      Eigen::VectorXd dPhi_ds_j = mGradientPhi_ds.row(j);
      auto dPhi_dx_j = dPhi_dr_j*drdx + dPhi_ds_j*dsdx;
      auto dPhi_dz_j = dPhi_dr_j*drdz + dPhi_ds_j*dsdz;
            
      mElementStiffnessMatrix(i,j) = detJ*mIntegrationWeights.dot((velocity.array().pow(2) * dPhi_dx_i.array() * dPhi_dx_j.array()).matrix()) +
        detJ*mIntegrationWeights.dot((velocity.array().pow(2) * dPhi_dz_i.array() * dPhi_dz_j.array()).matrix());
    }
  }
    
}

Eigen::MatrixXd AcousticTri::computeStiffnessTerm(const Eigen::MatrixXd &displacement) {

    Eigen::VectorXd Ku = mElementStiffnessMatrix*displacement.col(0);
    return Ku;
    
}

void AcousticTri::attachMaterialProperties(ExodusModel *model) {

    // Vp (m/s).
    mMaterialVelocityAtVertices = __attachMaterialProperties(model, "VP");

}

Eigen::MatrixXd AcousticTri::computeSourceTerm(double time) {

    // Initialize source vector (note: due to RVO I believe no memory re-allocation is occuring).
    Eigen::VectorXd F = Eigen::VectorXd::Zero(mNumIntPnt);
    
    return F;
}

void AcousticTri::computeSurfaceTerm() {

    std::cout << mElmNum << std::endl;

}

void AcousticTri::assembleElementMassMatrix(Mesh *mesh) {

    int i=0;
    Eigen::Matrix<double,2,2> Jinv;
    double detJ;
    Eigen::VectorXd elementMassMatrix(mNumIntPnt);
    std::tie(Jinv,detJ) = inverseJacobianAtPoint(0,0);
    elementMassMatrix = detJ*mIntegrationWeights;
    // assemble to shared nodes
    mesh->addFieldFromElement("m", mElmNum, mClsMap, elementMassMatrix);
    
}

void AcousticTri::setupTest(Mesh* mesh, Options options) {
    Eigen::VectorXd pts_x,pts_z;
    std::tie(pts_x,pts_z) = buildNodalPoints();
    // push nodal locations to shared dofs
    setInitialCondition(mesh,pts_x,pts_z,
                        options.IC_SquareSide_L(),
                        options.IC_Center_x(),options.IC_Center_z());
}

void AcousticTri::setInitialCondition(Mesh* mesh, Eigen::VectorXd& pts_x,Eigen::VectorXd& pts_z,
                                      double L, double x0, double z0) {
    
    double Lx = L;
    double Lz = L;
    const double PI = std::atan(1.0)*4.0;
    Eigen::VectorXd un = (PI/Lx*(pts_x.array()-(x0+L/2))).sin()*(PI/Lz*(pts_z.array()-(z0+L/2))).sin();
    Eigen::VectorXd vn = 0*pts_x;
    Eigen::VectorXd an = 0*pts_x;    
    mesh->setFieldFromElement("u", mElmNum, mClsMap, un);
    mesh->setFieldFromElement("v", mElmNum, mClsMap, vn);
    mesh->setFieldFromElement("a_", mElmNum, mClsMap, an);
    
}

double AcousticTri::checkTest(Mesh* mesh, Options options, const Eigen::MatrixXd &displacement, double time) {

    auto u_current = displacement.col(0);
    // exact solution
    Eigen::VectorXd pts_x,pts_z;
    std::tie(pts_x,pts_z) = buildNodalPoints();
    auto un_exact = exactSolution(pts_x,pts_z,
                                  options.IC_SquareSide_L(),
                                  options.IC_Center_x(),options.IC_Center_z(),
                                  time);

    auto element_error = (un_exact - u_current).array().abs().maxCoeff();

    return element_error;
}

Eigen::VectorXd AcousticTri::exactSolution(Eigen::VectorXd& pts_x,Eigen::VectorXd& pts_z,
                                           double L, double x0, double z0, double time) {

    const double PI = std::atan(1.0)*4.0;
    double Lx = L;
    double Lz = L;
    Eigen::VectorXd un_xz = (PI/Lx*(pts_x.array()-(x0+L/2))).sin()*(PI/Lz*(pts_z.array()-(z0+L/2))).sin();
    double velocity = mMaterialVelocityAtVertices.mean();
    double un_t = cos(PI/Lx*sqrt(2)*time*velocity);
    return un_t*un_xz;
}
