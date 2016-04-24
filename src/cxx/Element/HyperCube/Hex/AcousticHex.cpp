#include "AcousticHex.h"

AcousticHex::AcousticHex(Options options): Hexahedra(options) {

  // Allocate element vectors.
  mMssMat.setZero(mNumIntPnt);

  // Strain matrix.
  mElementStrain.setZero(3, mNumIntPnt);

}

Eigen::MatrixXd AcousticHex::computeStiffnessTerm(const Eigen::MatrixXd &displacement) {

  // Current gll point.
  int itr = 0;

  // Data structures we'll need here. Static arrays are allocated for free.
  // TODO: Look into a better way to deal with the temporarily allocated vectors. I believe that due to
  // TODO: RVO the integratedStiffnessMatrix should be ok, but perhaps jacobian_determinant could be handled better.
  Eigen::MatrixXd velocity_gradient(3, mNumIntPnt);
  Eigen::VectorXd detJ(mNumIntPnt);
  Eigen::VectorXd integratedStiffnessMatrix(mNumIntPnt);
  Eigen::Matrix<double,3,3> Jinv;
  double detJi;

  // Loop over all GLL points once to calculate the gradient of the pressure (u).
  Vector3d elementStrainRST(3);
  
  for (auto t_index = 0; t_index < mNumIntPtsT; t_index++) {
    for (auto s_index = 0; s_index < mNumIntPtsS; s_index++) {
      for (auto r_index = 0; r_index < mNumIntPtsR; r_index++) {

        // r,s,t coordinates.
        double t = mIntCoordT[t_index];
        double s = mIntCoordS[s_index];
        double r = mIntCoordR[r_index];

        // Get Jacobian determinant and its inverse
        std::tie(Jinv,detJi) = inverseJacobianAtPoint(r,s,t);
        detJ[itr] = detJi;

        // Calculate gradient in (r,s,t).
        elementStrainRST <<
          mGrd.row(r_index).dot(rVectorStride(displacement, s_index,t_index)),
          mGrd.row(s_index).dot(sVectorStride(displacement, r_index,t_index)),
          mGrd.row(t_index).dot(tVectorStride(displacement, r_index,s_index));
        
        // (r,s,t) -> (x,y,z) and save for kernel calculations.
        mElementStrain.col(itr) = Jinv*elementStrainRST;
        
        // interpolate material parameters at this node.
        double velocity = interpolateAtPoint(r,s,t).dot(mMaterialVelocityAtVertices);
        
        // apply material velocity (note v^2)
        velocity_gradient.col(itr) = velocity * velocity * mElementStrain.col(itr);
        
        itr++;

      }
    }
  }

  // Loop over all gll points again. Apply shape function derivative and integrate.
  itr = 0;
    
  for (auto t_index = 0; t_index < mNumIntPtsT; t_index++) {
    for (auto s_index = 0; s_index < mNumIntPtsS; s_index++) {
      for (auto r_index = 0; r_index < mNumIntPtsR; r_index++) {

        // r,s,t coordinates.
        double t = mIntCoordT[t_index];
        double s = mIntCoordS[s_index];
        double r = mIntCoordR[r_index];

        double detJi;
        std::tie(Jinv,detJi) = inverseJacobianAtPoint(r,s,t);


        // map reference gradient (lr,ls) to this element (lx,lz)
        auto lr = mGrd.col(r_index);
        auto ls = mGrd.col(s_index);
        auto lt = mGrd.col(t_index);
        
        auto dphi_r_dux = mIntWeightsS(s_index) * mIntWeightsT(t_index) *
          mIntWeightsR.dot(((rVectorStride(detJ, s_index, t_index)).array() *
                            rVectorStride(velocity_gradient.row(0), s_index,t_index).array() *
                            lr.array()).matrix());
        auto dphi_s_dux = mIntWeightsR(r_index) * mIntWeightsT(t_index) *
          mIntWeightsS.dot(((sVectorStride(detJ, r_index, t_index)).array() *
                            sVectorStride(velocity_gradient.row(0), r_index,t_index).array() *
                            ls.array()).matrix());
        auto dphi_t_dux = mIntWeightsR(r_index) * mIntWeightsS(s_index) *
          mIntWeightsT.dot(((tVectorStride(detJ, r_index, s_index)).array() *
                            tVectorStride(velocity_gradient.row(0), r_index,s_index).array() *
                            lt.array()).matrix());

        auto dphi_r_duy = mIntWeightsS(s_index) * mIntWeightsT(t_index) *
          mIntWeightsR.dot(((rVectorStride(detJ, s_index, t_index)).array() *
                            rVectorStride(velocity_gradient.row(1), s_index,t_index).array() *
                            lr.array()).matrix());
        auto dphi_s_duy = mIntWeightsR(r_index) * mIntWeightsT(t_index) *
          mIntWeightsS.dot(((sVectorStride(detJ, r_index, t_index)).array() *
                            sVectorStride(velocity_gradient.row(1), r_index,t_index).array() *
                            ls.array()).matrix());
        auto dphi_t_duy = mIntWeightsR(r_index) * mIntWeightsS(s_index) *
          mIntWeightsT.dot(((tVectorStride(detJ, r_index, s_index)).array() *
                            tVectorStride(velocity_gradient.row(1), r_index,s_index).array() *
                            lt.array()).matrix());
        
        auto dphi_r_duz = mIntWeightsS(s_index) * mIntWeightsT(t_index) *
          mIntWeightsR.dot(((rVectorStride(detJ, s_index, t_index)).array() *
                            rVectorStride(velocity_gradient.row(2), s_index,t_index).array() *
                            lr.array()).matrix());
        auto dphi_s_duz = mIntWeightsR(r_index) * mIntWeightsT(t_index) *
          mIntWeightsS.dot(((sVectorStride(detJ, r_index, t_index)).array() *
                            sVectorStride(velocity_gradient.row(2), r_index,t_index).array() *
                            ls.array()).matrix());
        auto dphi_t_duz = mIntWeightsR(r_index) * mIntWeightsS(s_index) *
          mIntWeightsT.dot(((tVectorStride(detJ, r_index, s_index)).array() *
                            tVectorStride(velocity_gradient.row(2), r_index,s_index).array() *
                            lt.array()).matrix());
        
        VectorXd dphi_rst_dux(3);
        dphi_rst_dux <<
          dphi_r_dux,
          dphi_s_dux,
          dphi_t_dux;

        VectorXd dphi_rst_duy(3);
        dphi_rst_duy <<
          dphi_r_duy,
          dphi_s_duy,
          dphi_t_duy;

        VectorXd dphi_rst_duz(3);
        dphi_rst_duz <<
          dphi_r_duz,
          dphi_s_duz,
          dphi_t_duz;
        
        integratedStiffnessMatrix(itr) =
          Jinv.row(0).dot(dphi_rst_dux) +
          Jinv.row(1).dot(dphi_rst_duy) +
          Jinv.row(2).dot(dphi_rst_duz);
        
        itr++;            
        
      }
            
    }
  }

  return integratedStiffnessMatrix;
}

void AcousticHex::attachMaterialProperties(ExodusModel *model) {

  // Vp (m/s).
  mMaterialVelocityAtVertices = __attachMaterialProperties(model, "VP");

}

Eigen::MatrixXd AcousticHex::computeSourceTerm(double time) {

  std::cerr << "ComputeSourceTerm: Not implemented yet!\n";
  exit(1);
  // // Initialize source vector (note: due to RVO I believe no memory re-allocation is occuring).
    // Eigen::VectorXd F = Eigen::VectorXd::Zero(mNumIntPnt);

    // // For all sources tagging along with this element.
    // for (auto &source: mSrc) {

    //     // TODO: May make this more efficient (i.e. allocation every loop)?
    //     // Evaluate shape functions at source (eps, eta). Save the lagrange coefficients in current_source.
    //     Eigen::VectorXd current_source = interpolateLagrangePolynomials(
    //             source->ReferenceLocationEps(), source->ReferenceLocationEta(), mPlyOrd);

    //     // Loop over gll points
    //     for (auto eta_index = 0; eta_index < mNumIntPtsEta; eta_index++) {
    //         for (auto eps_index = 0; eps_index < mNumIntPtsEps; eps_index++) {

    //             double eps = mIntCrdEps[eps_index];
    //             double eta = mIntCrdEta[eta_index];

    //             double detJ;
    //             Eigen::Matrix2d Jinv;
    //             std::tie(Jinv, detJ) = inverseJacobianAtPoint(eps,eta);

    //             // Calculate the coefficients needed to integrate to the delta function.
    //             current_source[eps_index + eta_index*mNumIntPtsEps] /=
    //                 (mIntWgtEps(eps_index) * mIntWgtEta(eta_index) * detJ);


    //         }
    //     }

    //     // Scale by the source amplitude.
    //     current_source *= source->fire(time);
    //     F += current_source;

    // }

    // return F;
}

void AcousticHex::computeSurfaceTerm() {

  std::cerr << "ComputeSurfaceTerm: Not implemented yet!\n";
  exit(1);

}

void AcousticHex::assembleElementMassMatrix(Mesh *mesh) {

  int i=0;
  Eigen::Matrix<double,3,3> Jinv;
  double detJ;
  Eigen::VectorXd elementMassMatrix(mNumIntPnt);
  for (auto t_index = 0; t_index < mNumIntPtsT; t_index++) {
    for (auto s_index = 0; s_index < mNumIntPtsS; s_index++) {
      for (auto r_index = 0; r_index < mNumIntPtsR; r_index++) {
        double r = mIntCoordR[r_index];
        double s = mIntCoordS[s_index];
        double t = mIntCoordT[t_index];
        std::tie(Jinv,detJ) = inverseJacobianAtPoint(r,s,t);
        elementMassMatrix[i] = detJ * mIntWeightsR[r_index] * mIntWeightsS[s_index] * mIntWeightsT[t_index];
        i++;
      }
    }
  }
  // assemble to shared nodes
  mesh->addFieldFromElement("m", mElmNum, mClsMap, elementMassMatrix);
}

double AcousticHex::checkTest(Mesh* mesh, Options options, const Eigen::MatrixXd &displacement, double time) {

  auto u_current = displacement.col(0);
  // exact solution
  Eigen::VectorXd pts_x,pts_y,pts_z;
  std::tie(pts_x,pts_y,pts_z) = buildNodalPoints();
  auto un_exact = exactSolution(pts_x,pts_y,pts_z,
                                options.IC_SquareSide_L(),
                                options.IC_Center_x(),options.IC_Center_y(),options.IC_Center_z(),
                                time);

  // mesh->setFieldFromElement("u_exact", mElmNum, mClsMap, un_exact);
  
  auto element_error = (un_exact - u_current).array().abs().maxCoeff();

  return element_error;
}

void AcousticHex::setupTest(Mesh* mesh, Options options) {
  Eigen::VectorXd pts_x,pts_y,pts_z;
  std::tie(pts_x,pts_y,pts_z) = buildNodalPoints();
  // push nodal locations to shared dofs
  setInitialCondition(mesh,pts_x,pts_y,pts_z,
                      options.IC_SquareSide_L(),
                      options.IC_Center_x(),options.IC_Center_y(),options.IC_Center_z());
}

void AcousticHex::setInitialCondition(Mesh* mesh,
                                       Eigen::VectorXd& pts_x,
                                       Eigen::VectorXd& pts_y,
                                       Eigen::VectorXd& pts_z,
                                       double L, double x0, double y0, double z0) {
    
  double Lx = L;
  double Ly = L;
  double Lz = L;
  const double PI = std::atan(1.0)*4.0;
  // exact solution
  Eigen::VectorXd un =
    (PI/Lx*(pts_x.array()-(x0+L/2))).sin()
    *(PI/Ly*(pts_y.array()-(y0+L/2))).sin()
    *(PI/Lz*(pts_z.array()-(z0+L/2))).sin();
  
  Eigen::VectorXd vn = 0*pts_x;
  Eigen::VectorXd an = 0*pts_x;    
  mesh->setFieldFromElement("u", mElmNum, mClsMap, un);
  mesh->setFieldFromElement("v", mElmNum, mClsMap, vn);
  mesh->setFieldFromElement("a_", mElmNum, mClsMap, an);
    
}

Eigen::VectorXd AcousticHex::exactSolution(Eigen::VectorXd& pts_x,
                                           Eigen::VectorXd& pts_y,
                                           Eigen::VectorXd& pts_z,
                                           double L, double x0, double y0, double z0, double time) {

  const double PI = std::atan(1.0)*4.0;
  double Lx = L;
  double Ly = L;
  double Lz = L;
  Eigen::VectorXd un_xyz =
    (PI/Lx*(pts_x.array()-(x0+L/2))).sin()
    *(PI/Ly*(pts_y.array()-(y0+L/2))).sin()
    *(PI/Lz*(pts_z.array()-(z0+L/2))).sin();  
  double velocity = mMaterialVelocityAtVertices.mean();
  double un_t = cos(PI/Lx*sqrt(3)*(time)*velocity);
  return un_t*un_xyz;
}

void AcousticHex::attachReceiver(std::vector<std::shared_ptr<Receiver>> &receivers) {
  printf("TODO: attachReciever");
  exit(1);
}

MatrixXd AcousticHex::interpolateFieldAtPoint(const VectorXd &pnt) {
  MatrixXd field;
  printf("TODO: interpolateFieldAtPoint");
  exit(1);
  return field;
}
void AcousticHex::recordField(const MatrixXd &u) {
  printf("TODO: recordField");
  exit(1);
}

