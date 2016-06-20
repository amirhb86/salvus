#include <Mesh/Mesh.h>
#include <Source/Source.h>
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>
#include <Receiver/Receiver.h>
#include <Element/HyperCube/Quad.h>
#include <Element/HyperCube/QuadP1.h>

extern "C" {
#include <Element/HyperCube/Autogen/quad_autogen.h>
}

using namespace Eigen;

template <typename ConcreteShape>
Quad<ConcreteShape>::Quad(std::unique_ptr<Options> const &options) {

  /* Ensure we've set parameters correctly. */
  assert(options->PolynomialOrder() > 0);

  mPlyOrd = options->PolynomialOrder();
  mNumDofVtx = 1;
  mNumDofEdg = mPlyOrd - 1;
  mNumDofFac = (mPlyOrd - 1) * (mPlyOrd - 1);
  mNumDofVol = 0;

  mGrd = Quad<ConcreteShape>::setupGradientOperator(mPlyOrd);
  mClsMap = Quad<ConcreteShape>::ClosureMappingForOrder(mPlyOrd);
  mIntCrdR = Quad<ConcreteShape>::GllPointsForOrder(mPlyOrd);
  mIntCrdS = Quad<ConcreteShape>::GllPointsForOrder(mPlyOrd);
  mIntWgtR = Quad<ConcreteShape>::GllIntegrationWeightsForOrder(mPlyOrd);
  mIntWgtS = Quad<ConcreteShape>::GllIntegrationWeightsForOrder(mPlyOrd);

  mNumIntPtsS = mIntCrdS.size();
  mNumIntPtsR = mIntWgtR.size();
  mNumIntPnt = mNumIntPtsS * mNumIntPtsR;

  mDetJac.setZero(mNumIntPnt);
  mParWork.setZero(mNumIntPnt);
  mStiffWork.setZero(mNumIntPnt);
  mGradWork.setZero(mNumIntPnt, mNumDim);

}

template <typename ConcreteShape>
VectorXd Quad<ConcreteShape>::GllPointsForOrder(const int order) {
  VectorXd gll_points(order + 1);
  if (order == 1) {
    gll_coordinates_order1_square(gll_points.data());
  } else if (order == 2) {
    gll_coordinates_order2_square(gll_points.data());
  } else if (order == 3) {
    gll_coordinates_order3_square(gll_points.data());
  } else if (order == 4) {
    gll_coordinates_order4_square(gll_points.data());
  } else if (order == 5) {
    gll_coordinates_order5_square(gll_points.data());
  } else if (order == 6) {
    gll_coordinates_order6_square(gll_points.data());
  } else if (order == 7) {
    gll_coordinates_order7_square(gll_points.data());
  } else if (order == 8) {
    gll_coordinates_order8_square(gll_points.data());
  } else if (order == 9) {
    gll_coordinates_order9_square(gll_points.data());
  } else if (order == 10) {
    gll_coordinates_order10_square(gll_points.data());
  }
  return gll_points;
}

template <typename ConcreteShape>
VectorXd Quad<ConcreteShape>::GllIntegrationWeightsForOrder(const int order) {
  VectorXd integration_weights(order + 1);
  if (order == 1) {
    gll_weights_order1_square(integration_weights.data());
  } else if (order == 2) {
    gll_weights_order2_square(integration_weights.data());
  } else if (order == 3) {
    gll_weights_order3_square(integration_weights.data());
  } else if (order == 4) {
    gll_weights_order4_square(integration_weights.data());
  } else if (order == 5) {
    gll_weights_order5_square(integration_weights.data());
  } else if (order == 6) {
    gll_weights_order6_square(integration_weights.data());
  } else if (order == 7) {
    gll_weights_order7_square(integration_weights.data());
  } else if (order == 8) {
    gll_weights_order8_square(integration_weights.data());
  } else if (order == 9) {
    gll_weights_order9_square(integration_weights.data());
  } else if (order == 10) {
    gll_weights_order10_square(integration_weights.data());
  }
  return integration_weights;
}

template <typename ConcreteShape>
VectorXi Quad<ConcreteShape>::ClosureMappingForOrder(const int order) {
  VectorXi closure_mapping((order + 1) * (order + 1));
  if (order == 1) {
    closure_mapping_order1_square(closure_mapping.data());
  } else if (order == 2) {
    closure_mapping_order2_square(closure_mapping.data());
  } else if (order == 3) {
    closure_mapping_order3_square(closure_mapping.data());
  } else if (order == 4) {
    closure_mapping_order4_square(closure_mapping.data());
  } else if (order == 5) {
    closure_mapping_order5_square(closure_mapping.data());
  } else if (order == 6) {
    closure_mapping_order6_square(closure_mapping.data());
  } else if (order == 7) {
    closure_mapping_order7_square(closure_mapping.data());
  } else if (order == 8) {
    closure_mapping_order8_square(closure_mapping.data());
  } else if (order == 9) {
    closure_mapping_order9_square(closure_mapping.data());
  } else if (order == 10) {
    closure_mapping_order10_square(closure_mapping.data());
  }
  return closure_mapping;
}

template <typename ConcreteShape>
VectorXd Quad<ConcreteShape>::rVectorStride(const Ref<const VectorXd>& f, const int s_ind,
                                               const int numPtsS, const int numPtsR) {
  return Map<const VectorXd> (f.data() + s_ind * numPtsS, numPtsR);
}

template <typename ConcreteShape>
VectorXd Quad<ConcreteShape>::sVectorStride(const Ref<const VectorXd>& f, const int r_ind,
                                               const int numPtsS, const int numPtsR) {
  return Map<const VectorXd, 0, InnerStride<>> (f.data() + r_ind, numPtsS, InnerStride<>(numPtsR));
}


template <typename ConcreteShape>
void Quad<ConcreteShape>::attachVertexCoordinates(std::unique_ptr<Mesh> const &mesh) {

  Vec coordinates_local;      /* TODO:: DEALLOCATE THESE? */
  PetscInt coordinate_buffer_size;
  PetscSection coordinate_section;
  PetscReal *coordinates_buffer = NULL;

  DMGetCoordinatesLocal(mesh->DistributedMesh(), &coordinates_local);
  DMGetCoordinateSection(mesh->DistributedMesh(), &coordinate_section);
  DMPlexVecGetClosure(mesh->DistributedMesh(), coordinate_section, coordinates_local, mElmNum,
                      &coordinate_buffer_size, &coordinates_buffer);
  std::vector<PetscReal> coordinates_element(coordinates_buffer, coordinates_buffer + coordinate_buffer_size);
  DMPlexVecRestoreClosure(mesh->DistributedMesh(), coordinate_section, coordinates_local, mElmNum,
                          &coordinate_buffer_size, &coordinates_buffer);

  // Get all coordinates.
  for (int i = 0; i < mNumVtx; i++) {
    mVtxCrd(i,0) = coordinates_element[mNumDim * i + 0];
    mVtxCrd(i,1) = coordinates_element[mNumDim * i + 1];
  }

  // Save edge maps.
  mEdgMap = mesh->EdgeNumbers(mElmNum);

  // Save element center
  mElmCtr << mVtxCrd.col(0).mean(),
             mVtxCrd.col(1).mean();

}

template <typename ConcreteShape>
double Quad<ConcreteShape>::CFL_constant() {
  if(mPlyOrd == 3) {
    return 3.0; // determined by hand (about 10% conservative)
  }
  if(mPlyOrd == 4) {
    return 2.0; // ??determined by hand (about 10% conservative)
  }
  else {
    std::cerr << "ERROR: Order CFL_constant not implemented yet\n";
    exit(1);
  }
}

template <typename ConcreteShape>
double Quad<ConcreteShape>::estimatedElementRadius() {

  Matrix2d invJ;
  double detJ;
  
  Matrix2d invJac;
  Vector2d refGrad;
  int num_pts = mNumIntPtsR*mNumIntPtsS;
  VectorXd h_pts(num_pts);
  
  // Loop over all GLL points.
  for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
    for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

      // gll index.
      int index = r_ind + s_ind * mNumIntPtsR;

      // (r,s,t) coordinates for this point.
      double r = mIntCrdR(r_ind);
      double s = mIntCrdS(s_ind);

      // Optimized gradient for tensorized GLL basis.
      std::tie(invJ, detJ) = ConcreteShape::inverseJacobianAtPoint(r, s, mVtxCrd);
      Matrix2d J = invJ.inverse();
      VectorXcd eivals = J.eigenvalues();
      // get minimum h (smallest direction)
      Vector2d eivals_norm;
      for(int i=0;i<2;i++) {
        eivals_norm(i) = std::norm(eivals[i]);
      }
      h_pts(index) = eivals_norm.minCoeff();
    
    }
  }
  return h_pts.minCoeff();
  
}


template <typename ConcreteShape>
void Quad<ConcreteShape>::attachMaterialProperties(std::unique_ptr<ExodusModel> const &model, std::string parameter) {
  Vector4d material_at_vertices;
  for (int i = 0; i < mNumVtx; i++) {
    material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(
        mElmCtr, parameter, i);
  }
  mPar[parameter] = material_at_vertices;
}

template <typename ConcreteShape>
void Quad<ConcreteShape>::attachReceiver(std::vector<std::unique_ptr<Receiver>> receivers) {

  for (auto &rec: receivers) {
    double x1 = rec->PysLocX1();
    double x2 = rec->PysLocX2();
    if (ConcreteShape::checkHull(x1, x2, mVtxCrd)) {
      Vector2d ref_loc = ConcreteShape::inverseCoordinateTransform(x1, x2, mVtxCrd);
      rec->SetRefLocR(ref_loc(0));
      rec->SetRefLocS(ref_loc(1));
      mRec.push_back(std::move(rec));
    }
  }
}

template <typename ConcreteShape>
bool Quad<ConcreteShape>::attachSource(std::unique_ptr<Source> &source, const bool finalize) {
  if (!source) { return false; }
  double x1 = source->LocX();
  double x2 = source->LocZ();
  if (ConcreteShape::checkHull(x1, x2, mVtxCrd)) {
    if (!finalize) { return true; }
    Vector2d ref_loc = ConcreteShape::inverseCoordinateTransform(x1, x2, mVtxCrd);
    source->SetLocR(ref_loc(0));
    source->SetLocS(ref_loc(1));
    mSrc.push_back(std::move(source));
    return true;
  }
  return false;
}

template <typename ConcreteShape>
VectorXd Quad<ConcreteShape>::getDeltaFunctionCoefficients(const double r, const double s) {

  Matrix2d _;
  mParWork = interpolateLagrangePolynomials(r, s, mPlyOrd);
  for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
    for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

      double ri = mIntCrdR(r_ind);
      double si = mIntCrdS(s_ind);

      double detJac;
      std::tie(_, detJac) = ConcreteShape::inverseJacobianAtPoint(ri, si, mVtxCrd);

      mParWork(r_ind + s_ind * mNumIntPtsR) /=
          (mIntWgtR(r_ind) * mIntWgtS(s_ind) * detJac);

    }
  }
  return mParWork;
}

template <typename ConcreteShape>
VectorXd Quad<ConcreteShape>::interpolateLagrangePolynomials(const double r, const double s,
                                                          const int order) {

  assert(order > 0 && order < 11);

  int n_points = (order + 1) * (order + 1);
  VectorXd gll_coeffs(n_points);
  if (order == 1) {
    interpolate_order1_square(r, s, gll_coeffs.data());
  } else if (order == 2) {
    interpolate_order2_square(r, s, gll_coeffs.data());
  } else if (order == 3) {
    interpolate_order3_square(r, s, gll_coeffs.data());
  } else if (order == 4) {
    interpolate_order4_square(r, s, gll_coeffs.data());
  } else if (order == 5) {
    interpolate_order5_square(r, s, gll_coeffs.data());
  } else if (order == 6) {
    interpolate_order6_square(r, s, gll_coeffs.data());
  } else if (order == 7) {
    interpolate_order7_square(r, s, gll_coeffs.data());
  } else if (order == 8) {
    interpolate_order8_square(r, s, gll_coeffs.data());
  } else if (order == 9) {
    interpolate_order9_square(r, s, gll_coeffs.data());
  } else if (order == 10) {
    interpolate_order10_square(r, s, gll_coeffs.data());
  }
  return gll_coeffs;
}


template <typename ConcreteShape>
MatrixXd Quad<ConcreteShape>::setupGradientOperator(const int order) {

  int num_pts_r = Quad<ConcreteShape>::GllPointsForOrder(order).size();
  int num_pts_s = Quad<ConcreteShape>::GllPointsForOrder(order).size();
  double eta = Quad<ConcreteShape>::GllPointsForOrder(order)(0);

  MatrixXd grad(num_pts_s, num_pts_r);
  MatrixXd test(num_pts_s, num_pts_r);
  for (int i = 0; i < num_pts_r; i++) {
    double eps = Quad<ConcreteShape>::GllPointsForOrder(order)(i);
    if (order == 1) {
      interpolate_eps_derivative_order1_square(eta, test.data());
    } else if (order == 2) {
      interpolate_eps_derivative_order2_square(eps, eta, test.data());
    } else if (order == 3) {
      interpolate_eps_derivative_order3_square(eps, eta, test.data());
    } else if (order == 4) {
      interpolate_eps_derivative_order4_square(eps, eta, test.data());
    } else if (order == 5) {
      interpolate_eps_derivative_order5_square(eps, eta, test.data());
    } else if (order == 6) {
      interpolate_eps_derivative_order6_square(eps, eta, test.data());
    } else if (order == 7) {
      interpolate_eps_derivative_order7_square(eps, eta, test.data());
    } else if (order == 8) {
      interpolate_eps_derivative_order8_square(eps, eta, test.data());
    } else if (order == 9) {
      interpolate_eps_derivative_order9_square(eps, eta, test.data());
    } else if (order == 10) {
      interpolate_eps_derivative_order10_square(eps, eta, test.data());
    }
    grad.row(i) = test.col(0);
  }
  return grad;
}

template <typename ConcreteShape>
MatrixXd Quad<ConcreteShape>::computeGradient(const Ref<const VectorXd> &field) {

  Matrix2d invJac;
  Vector2d refGrad;

  // Loop over all GLL points.
  for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
    for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

      // gll index.
      int index = r_ind + s_ind * mNumIntPtsR;

      // (r,s) coordinates for this point.
      double r = mIntCrdR(r_ind);
      double s = mIntCrdS(s_ind);
      std::tie(invJac, mDetJac(index)) = ConcreteShape::inverseJacobianAtPoint(r, s, mVtxCrd);
      
      // optimized version (17us->13us)
      refGrad.setZero(2);
      for(int i=0;i<mNumIntPtsR;i++) {
        refGrad(0) += mGrd(r_ind,i)*field(i + s_ind * mNumIntPtsR);
        refGrad(1) += mGrd(s_ind,i)*field(r_ind + i * mNumIntPtsR);
      }
      mGradWork.row(index).noalias() = invJac * refGrad;        
      
    }
  }

  return mGradWork;

}

template <typename ConcreteShape>
VectorXd Quad<ConcreteShape>::ParAtIntPts(const std::string &par) {

  for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
    for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

      double r = mIntCrdR(r_ind);
      double s = mIntCrdS(s_ind);
      mParWork(r_ind + s_ind*mNumIntPtsR) = ConcreteShape::interpolateAtPoint(r,s).dot(mPar[par]);

    }
  }

  return mParWork;
}

template <typename ConcreteShape>
VectorXd Quad<ConcreteShape>::applyTestAndIntegrate(const Ref<const VectorXd> &f) {

  int i = 0;
  double detJac;
  Matrix2d invJac;
  VectorXd result(mNumIntPnt);
  for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
    for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

      // gll index.
      int index = r_ind + s_ind * mNumIntPtsR;

      // (r,s) coordinate at this point.
      double r = mIntCrdR(r_ind);
      double s = mIntCrdS(s_ind);

      std::tie(invJac,detJac) = ConcreteShape::inverseJacobianAtPoint(r,s,mVtxCrd);
      result(index) = f(index) * detJac * mIntWgtR(r_ind) * mIntWgtS(s_ind);

    }
  }

  return result;

}

template <typename ConcreteShape>
Eigen::Vector2d Quad<ConcreteShape>::getEdgeNormal(const PetscInt edg) {

  Vector2d n;
  double x0, x1, y0, y1;
  for (int i = 0; i < mNumVtx; i++) {
    if (mEdgMap[i] == edg) {
      x0 = mVtxCrd(i, 0);
      x1 = mVtxCrd((i+1) % mNumVtx, 0);
      y0 = mVtxCrd(i, 1);
      y1 = mVtxCrd((i+1) % mNumVtx, 1);
    }
  }

  n(0) = -1 * (y1 - y0);
  n(1) = +1 * (x1 - x0);
  return n / n.norm();

}

template <typename ConcreteShape>
Eigen::VectorXd Quad<ConcreteShape>::applyTestAndIntegrateEdge(const Eigen::Ref<const Eigen::VectorXd> &f,
                                                               const PetscInt edg) {

  // Allocate return vector (return-value optimization)?
  Eigen::VectorXd result = Eigen::VectorXd::Zero(mNumIntPnt);

  // get edge vertices.
  int start, stride;
  double x0, x1, y0, y1;
  for (int i = 0; i < mNumVtx; i++) {
    if (mEdgMap[i] == edg) {
      x0 = mVtxCrd(i, 0);
      x1 = mVtxCrd((i+1) % mNumVtx, 0);
      y0 = mVtxCrd(i, 1);
      y1 = mVtxCrd((i+1) % mNumVtx, 1);
      if      (i == 0) { start = 0; stride = 1; }
      else if (i == 1) { start = mNumIntPtsR-1; stride = mNumIntPtsR; }
      else if (i == 2) { start = mNumIntPtsR*mNumIntPtsS - 1; stride = -1; }
      else             { start = mNumIntPtsR * (mNumIntPtsS - 1); stride = -1 * mNumIntPtsR; }
    }
  }

  // compute 'edge jacobian'.
  double dx = x1 - x0;
  double dy = y1 - y0;
  double d = sqrt(dx*dx + dy*dy) / 2.0;

  // compute coefficients.
  for (int i = start, j = 0; j < mNumIntPtsR; i+=stride, j++) {
    result(i) = f(i) * d * mIntWgtR(j);
  }

  return result;

}

template <typename ConcreteShape>
VectorXd Quad<ConcreteShape>::applyGradTestAndIntegrate(const Ref<const MatrixXd> &f) {

  Matrix2d invJac;
  for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
    
    for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

      double r = mIntCrdR(r_ind);
      double s = mIntCrdS(s_ind);

      double _;
      std::tie(invJac, _) = ConcreteShape::inverseJacobianAtPoint(r,s,mVtxCrd);

      // double dphi_r_dfx = mIntWgtS(s_ind) *
      //     mIntWgtR.dot(((rVectorStride(mDetJac, s_ind, mNumIntPtsS, mNumIntPtsR)).array() *
      //     rVectorStride(f.col(0), s_ind, mNumIntPtsS, mNumIntPtsR).array() *
      //     mGrd.col(r_ind).array()).matrix());

      // double dphi_s_dfx = mIntWgtR(r_ind) *
      //     mIntWgtS.dot(((sVectorStride(mDetJac, r_ind, mNumIntPtsS, mNumIntPtsR)).array() *
      //     sVectorStride(f.col(0), r_ind, mNumIntPtsS, mNumIntPtsR).array() *
      //     mGrd.col(s_ind).array()).matrix());

      // double dphi_r_dfy = mIntWgtS(s_ind) *
      //     mIntWgtR.dot(((rVectorStride(mDetJac, s_ind, mNumIntPtsS, mNumIntPtsR)).array() *
      //     rVectorStride(f.col(1), s_ind, mNumIntPtsS, mNumIntPtsR).array() *
      //     mGrd.col(r_ind).array()).matrix());

      // double dphi_s_dfy = mIntWgtR(r_ind) *
      //     mIntWgtS.dot(((sVectorStride(mDetJac, r_ind, mNumIntPtsS, mNumIntPtsR)).array() *
      //     sVectorStride(f.col(1), r_ind, mNumIntPtsS, mNumIntPtsR).array() *
      //     mGrd.col(s_ind).array()).matrix());
      
      auto lr = mGrd.col(r_ind);
      auto ls = mGrd.col(s_ind);
      
      double dphi_r_dfx = 0;
      double dphi_s_dfx = 0;
      
      double dphi_r_dfy = 0;
      double dphi_s_dfy = 0;
      
      for(int i=0;i<mNumIntPtsR;i++) {
        int r_index = i + s_ind * mNumIntPtsR;
        int s_index = r_ind + i * mNumIntPtsR;

        dphi_r_dfx += mDetJac[r_index] * f(r_index,0) * lr[i] * mIntWgtR[i];
        dphi_s_dfx += mDetJac[s_index] * f(s_index,0) * ls[i] * mIntWgtS[i];

        // -- y --
        dphi_r_dfy += mDetJac[r_index] * f(r_index,1) * lr[i] * mIntWgtR[i];
        dphi_s_dfy += mDetJac[s_index] * f(s_index,1) * ls[i] * mIntWgtS[i];

      }

      dphi_r_dfx *= mIntWgtS(s_ind);
      dphi_s_dfx *= mIntWgtR(r_ind);
      
      dphi_r_dfy *= mIntWgtS(s_ind);
      dphi_s_dfy *= mIntWgtR(r_ind);
      
      Vector2d dphi_epseta_dfx, dphi_epseta_dfy;
      dphi_epseta_dfx << dphi_r_dfx, dphi_s_dfx;
      dphi_epseta_dfy << dphi_r_dfy, dphi_s_dfy;

      mStiffWork(r_ind + s_ind * mNumIntPtsR) =
        invJac.row(0).dot(dphi_epseta_dfx) +
        invJac.row(1).dot(dphi_epseta_dfy);

    }
  }

  return mStiffWork;

}

template <typename ConcreteShape>
double Quad<ConcreteShape>::integrateField(const Eigen::Ref<const Eigen::VectorXd> &field) {

  double val = 0;
  Matrix2d inverse_Jacobian;
  double detJ;
  for (int i = 0; i < mNumIntPtsS; i++) {
    for (int j = 0; j < mNumIntPtsR; j++) {

      double r = mIntCrdR(j);
      double s = mIntCrdS(i);
      std::tie(inverse_Jacobian, detJ) = ConcreteShape::inverseJacobianAtPoint(r, s, mVtxCrd);
      val += field(j + i * mNumIntPtsR) * mIntWgtR(j) *
          mIntWgtS(i) * detJ;

    }
  }

  return val;
}

template <typename ConcreteShape>
void Quad<ConcreteShape>::setBoundaryConditions(std::unique_ptr<Mesh> const &mesh) {

  mBndElm = false;
  for (auto &keys: mesh ->BoundaryElementFaces()) {
    auto boundary_name = keys.first;
    auto element_in_boundary = keys.second;
    if (element_in_boundary.find(mElmNum) != element_in_boundary.end()) {
      mBndElm = true;
      mBnd[boundary_name] = element_in_boundary[mElmNum];
    }
  }

}

template <typename ConcreteShape>
void Quad<ConcreteShape>::applyDirichletBoundaries(std::unique_ptr<Mesh> const &mesh,
                                                   std::unique_ptr<Options> const &options,
                                                   const std::string &fieldname) {

  if (! mBndElm) return;

  double value = 0;
  auto dirchlet_boundary_names = options->DirichletBoundaries();
  for (auto &bndry: dirchlet_boundary_names) {
    auto faceids = mBnd[bndry];
    for (auto &faceid: faceids) {
      auto field = mesh->getFieldOnFace(fieldname, faceid);
      field = 0 * field.array() + value;
      mesh->setFieldFromFace(fieldname, faceid, field);
    }
  }
}

// Instantiate combinatorical cases.
template class Quad<QuadP1>;

