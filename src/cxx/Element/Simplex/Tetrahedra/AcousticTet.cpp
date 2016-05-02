#include "AcousticTet.h"

AcousticTet::AcousticTet(Options options): Tetrahedra(options) {
    
  // Allocate element vectors.
  mMssMat.setZero(mNumIntPnt);

  // Strain matrix.
  mElementStrain.setZero(3, mNumIntPnt);

}

void AcousticTet::buildStiffnessMatrix() {

  Eigen::Matrix3d invJ;
  double detJ;
  std::tie(invJ,detJ) = inverseJacobianAtPoint(0,0,0);
  //Jinv= rx, sx, tx,
  //      ry, sy, ty,
  //      rz, sz, tz;
  auto drdx = invJ(0,0);
  auto dsdx = invJ(0,1);
  auto dtdx = invJ(0,2);

  auto drdy = invJ(1,0);
  auto dsdy = invJ(1,1);
  auto dtdy = invJ(1,2);

  auto drdz = invJ(2,0);
  auto dsdz = invJ(2,1);
  auto dtdz = invJ(2,2);
  
  // build material on all nodes
  Eigen::VectorXd velocity(mNumIntPnt);
  mElementStiffnessMatrix.resize(mNumIntPnt,mNumIntPnt);
    
  // interpolate velocity at all nodes
  for(int i=0;i<mNumIntPnt;i++) {
    auto r = mIntegrationCoordinates_r[i];
    auto s = mIntegrationCoordinates_s[i];
    auto t = mIntegrationCoordinates_t[i];
    velocity(i) = interpolateAtPoint(r,s,t).dot(mMaterialVelocityAtVertices);
  }
  // if(mElmNum == 0) {
  //   std::cout << "Velocity(k=0)=" << velocity.transpose() << "\n";
  // }
  MatrixXd Kxx(mNumIntPnt,mNumIntPnt);
  // loop over matrix(i,j)
  for(int i=0;i<mNumIntPnt;i++) {
      
    Eigen::VectorXd dPhi_dr_i = mGradientPhi_dr.row(i);
    Eigen::VectorXd dPhi_ds_i = mGradientPhi_ds.row(i);
    Eigen::VectorXd dPhi_dt_i = mGradientPhi_dt.row(i);
    auto dPhi_dx_i = dPhi_dr_i*drdx + dPhi_ds_i*dsdx + dPhi_dt_i*dtdx;
    auto dPhi_dy_i = dPhi_dr_i*drdy + dPhi_ds_i*dsdy + dPhi_dt_i*dtdy;
    auto dPhi_dz_i = dPhi_dr_i*drdz + dPhi_ds_i*dsdz + dPhi_dt_i*dtdz;
    for(int j=0;j<mNumIntPnt;j++) {
      Eigen::VectorXd dPhi_dr_j = mGradientPhi_dr.row(j);
      Eigen::VectorXd dPhi_ds_j = mGradientPhi_ds.row(j);
      Eigen::VectorXd dPhi_dt_j = mGradientPhi_dt.row(j);
      auto dPhi_dx_j = dPhi_dr_j*drdx + dPhi_ds_j*dsdx + dPhi_dt_j*dtdx;
      auto dPhi_dy_j = dPhi_dr_j*drdy + dPhi_ds_j*dsdy + dPhi_dt_j*dtdy;
      auto dPhi_dz_j = dPhi_dr_j*drdz + dPhi_ds_j*dsdz + dPhi_dt_j*dtdz;
      
      Kxx(i,j) = detJ*mIntegrationWeights.dot((dPhi_dx_i.array() * dPhi_dx_j.array()).matrix());
      mElementStiffnessMatrix(i,j) =
        // with velocity according to model
        detJ*mIntegrationWeights.dot((velocity.array().pow(2) * dPhi_dx_i.array() * dPhi_dx_j.array()).matrix()) +
        detJ*mIntegrationWeights.dot((velocity.array().pow(2) * dPhi_dy_i.array() * dPhi_dy_j.array()).matrix()) +
        detJ*mIntegrationWeights.dot((velocity.array().pow(2) * dPhi_dz_i.array() * dPhi_dz_j.array()).matrix());
      // with velocity=1.0
      // detJ*mIntegrationWeights.dot((dPhi_dx_i.array() * dPhi_dx_j.array()).matrix()) +
        // detJ*mIntegrationWeights.dot((dPhi_dy_i.array() * dPhi_dy_j.array()).matrix()) +
        // detJ*mIntegrationWeights.dot((dPhi_dz_i.array() * dPhi_dz_j.array()).matrix());
    }
  }

  // if(mElmNum == 0) {
  //   std::cout << "element 1\n";
  //   auto v1x = mVtxCrd(0,0);
  //   auto v2x = mVtxCrd(0,1);
  //   auto v3x = mVtxCrd(0,2);
  //   auto v4x = mVtxCrd(0,3);
  //   printf("avg x: %f\n",(v1x+v2x+v3x+v4x)/4.0);
  
  //   auto v1y = mVtxCrd(1,0);
  //   auto v2y = mVtxCrd(1,1);
  //   auto v3y = mVtxCrd(1,2);
  //   auto v4y = mVtxCrd(1,3);
  
  //   auto v1z = mVtxCrd(2,0);
  //   auto v2z = mVtxCrd(2,1);
  //   auto v3z = mVtxCrd(2,2);
  //   auto v4z = mVtxCrd(2,3);
  //   printf("v1=(%f,%f,%f),v2=(%f,%f,%f),v3=(%f,%f,%f),v4=(%f,%f,%f)\n",
  //          v1x,v1y,v1z,
  //          v2x,v2y,v2z,
  //          v3x,v3y,v3z,
  //          v4x,v4y,v4z);
  //   printf("center=(%f,%f,%f)\n",mElmCtr(0),mElmCtr(1),mElmCtr(2));
  //   std::cout << "invJ=" << invJ << "\n";
  //   std::cout << "|J|=" << detJ << "\n";
  //   std::cout << "K[0,:]=" << mElementStiffnessMatrix.row(0) << "\n";
  //   std::cout << "Kxx[0,:]=" << Kxx.row(0) << "\n";
  //   std::cout << "dPhi_dr.row(0)=" << mGradientPhi_dr.row(0) << "\n";
  //   std::cout << "dPhi_ds.row(0)=" << mGradientPhi_ds.row(0) << "\n";
  //   std::cout << "dPhi_dt.row(0)=" << mGradientPhi_dt.row(0) << "\n";
  //   std::cout << "wn=" << mIntegrationWeights.transpose() << "\n";
  //   // exit(1);
  // }
  
}

Eigen::MatrixXd AcousticTet::computeStiffnessTerm(const Eigen::MatrixXd &displacement) {

  Eigen::VectorXd Ku = mElementStiffnessMatrix*displacement.col(0);
  // if(mElmNum == 0) {
  //   std::cout << "un_k=\n" << displacement.col(0).transpose() << "\n";
  //   std::cout << "an_k=\n" << Ku.transpose() << "\n";
  // }
  return Ku;
    
}

void AcousticTet::attachMaterialProperties(ExodusModel *model) {

  // Vp (m/s).
  mMaterialVelocityAtVertices = __attachMaterialProperties(model, "VP");

}

Eigen::MatrixXd AcousticTet::computeSourceTerm(double time) {

  // Initialize source vector (note: due to RVO I believe no memory re-allocation is occuring).
  Eigen::VectorXd F = Eigen::VectorXd::Zero(mNumIntPnt);

  printf("TODO: computeSourceTerm\n");
  exit(1);
  
  return F;
}

void AcousticTet::computeSurfaceTerm() {

  printf("TODO: computeSurfaceTerm\n");
  exit(1);
  std::cout << mElmNum << std::endl;

}

void AcousticTet::assembleElementMassMatrix(Mesh *mesh) {

  int i=0;
  Eigen::Matrix<double,3,3> Jinv;
  double detJ;
  Eigen::VectorXd elementMassMatrix(mNumIntPnt);

  std::tie(Jinv,detJ) = inverseJacobianAtPoint(0,0,0);
  elementMassMatrix = detJ*mIntegrationWeights;
  if(mElmNum == 0) {
    std::cout << "m(k==0)=" << elementMassMatrix.transpose() << "\n";
  }
  // assemble to shared nodes  
  mesh->addFieldFromElement("m", mElmNum, mClsMap, elementMassMatrix);
    
}

void AcousticTet::setupTest(Mesh* mesh, Options options) {
  Eigen::VectorXd pts_x,pts_y,pts_z;
  std::tie(pts_x,pts_y,pts_z) = buildNodalPoints();
  // push nodal locations to shared dofs
  setInitialCondition(mesh,pts_x,pts_y,pts_z,
                      options.IC_SquareSide_L(),
                      options.IC_Center_x(),options.IC_Center_y(),options.IC_Center_z());
}

void AcousticTet::setInitialCondition(Mesh* mesh,
                                      Eigen::VectorXd& pts_x,
                                      Eigen::VectorXd& pts_y,
                                      Eigen::VectorXd& pts_z,
                                      double L, double x0, double y0, double z0) {
    
  double Lx = L;
  double Ly = L;
  double Lz = L;
  const double PI = std::atan(1.0)*4.0;
  Eigen::VectorXd un =
    (PI/Lx*(pts_x.array()-(x0+L/2))).sin() *
    (PI/Ly*(pts_y.array()-(y0+L/2))).sin() *
    (PI/Lz*(pts_z.array()-(z0+L/2))).sin();
  Eigen::VectorXd vn = 0*pts_x;
  Eigen::VectorXd an = 0*pts_x;
  // un = pts_x.array().pow(2) + pts_z.array().pow(2) + pts_y.array().pow(2);
  // un = pts_x.array().pow(1) + pts_z.array().pow(1) + pts_y.array().pow(1);
  // un = (PI/Lx*(pts_x.array()-(x0+L/2))).sin();
  // un = (PI/Lx*(pts_y.array()-(y0+L/2))).sin();
  // un = un.array()*(PI/Lx*(pts_z.array()-(z0+L/2))).sin();
  mesh->setFieldFromElement("u", mElmNum, mClsMap, un);
  mesh->setFieldFromElement("v", mElmNum, mClsMap, vn);
  mesh->setFieldFromElement("a_", mElmNum, mClsMap, an);
    
}

double AcousticTet::checkTest(Mesh* mesh, Options options, const Eigen::MatrixXd &displacement, double time) {

  auto u_current = displacement.col(0);
  // exact solution
  Eigen::VectorXd pts_x,pts_y,pts_z;
  std::tie(pts_x,pts_y,pts_z) = buildNodalPoints();
  auto un_exact = exactSolution(pts_x,pts_y,pts_z,
                                options.IC_SquareSide_L(),
                                options.IC_Center_x(),options.IC_Center_x(),options.IC_Center_z(),
                                time);

  auto element_error = (un_exact - u_current).array().abs().maxCoeff();

  return element_error;
}

Eigen::VectorXd AcousticTet::exactSolution(Eigen::VectorXd& pts_x,
                                           Eigen::VectorXd& pts_y,
                                           Eigen::VectorXd& pts_z,
                                           double L, double x0, double y0, double z0, double time) {

  const double PI = std::atan(1.0)*4.0;
  double Lx = L;
  double Ly = L;
  double Lz = L;
  Eigen::VectorXd un_xyz =
    (PI/Lx*(pts_x.array()-(x0+L/2))).sin() *
    (PI/Ly*(pts_y.array()-(y0+L/2))).sin() *
    (PI/Lz*(pts_z.array()-(z0+L/2))).sin();
  double velocity = mMaterialVelocityAtVertices.mean();
  double un_t = cos(PI/Lx*sqrt(3)*time*velocity);
  return un_t*un_xyz;
}
