#include <Element/HyperCube/Quad/QuadP1.h>
#include "QuadNew.h"

using namespace Eigen;

template <typename Derived>
QuadNew<Derived>::QuadNew(Options options) {

  mPlyOrd = options.PolynomialOrder();
  mNumDofVtx = 1;
  mNumDofEdg = mPlyOrd - 1;
  mNumDofFac = (mPlyOrd - 1) * (mPlyOrd - 1);
  mNumDofVol = 0;

  mGrd = QuadNew<Derived>::setupGradientOperator(mPlyOrd);
  mClsMap = QuadNew<Derived>::ClosureMappingForOrder(mPlyOrd);
  mIntCrdR = QuadNew<Derived>::GllPointsForOrder(mPlyOrd);
  mIntCrdS = QuadNew<Derived>::GllPointsForOrder(mPlyOrd);
  mIntWgtR = QuadNew<Derived>::GllIntegrationWeightsForOrder(mPlyOrd);
  mIntWgtS = QuadNew<Derived>::GllIntegrationWeightsForOrder(mPlyOrd);

  mNumIntPtsS = mIntCrdS.size();
  mNumIntPtsR = mIntWgtR.size();
  mNumIntPnt = mNumIntPtsS * mNumIntPtsR;

  mDetJac.setZero(mNumIntPnt);
  mParWork.setZero(mNumIntPnt);
  mStiffWork.setZero(mNumIntPnt);
  mGradWork.setZero(mNumIntPnt, mNumDim);

}

template <typename Derived>
VectorXd QuadNew<Derived>::GllPointsForOrder(const int order) {
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

template <typename Derived>
VectorXd QuadNew<Derived>::GllIntegrationWeightsForOrder(const int order) {
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

template <typename Derived>
VectorXi QuadNew<Derived>::ClosureMappingForOrder(const int order) {
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
VectorXd QuadNew<ConcreteShape>::rVectorStride(const Ref<const VectorXd>& f, const int s_ind,
                                               const int numPtsS, const int numPtsR) {
  return Map<const VectorXd> (f.data() + s_ind * numPtsS, numPtsR);
}

template <typename ConcreteShape>
VectorXd QuadNew<ConcreteShape>::sVectorStride(const Ref<const VectorXd>& f, const int r_ind,
                                               const int numPtsS, const int numPtsR) {
  return Map<const VectorXd, 0, InnerStride<>> (f.data() + r_ind, numPtsS, InnerStride<>(numPtsR));
}


template <typename Derived>
void QuadNew<Derived>::attachVertexCoordinates(DM &distributed_mesh) {

  Vec coordinates_local;
  PetscInt coordinate_buffer_size;
  PetscSection coordinate_section;
  PetscReal *coordinates_buffer = NULL;

  DMGetCoordinatesLocal(distributed_mesh, &coordinates_local);
  DMGetCoordinateSection(distributed_mesh, &coordinate_section);
  DMPlexVecGetClosure(distributed_mesh, coordinate_section, coordinates_local, mElmNum,
                      &coordinate_buffer_size, &coordinates_buffer);
  std::vector<PetscReal> coordinates_element(coordinates_buffer, coordinates_buffer + coordinate_buffer_size);
  DMPlexVecRestoreClosure(distributed_mesh, coordinate_section, coordinates_local, mElmNum,
                          &coordinate_buffer_size, &coordinates_buffer);

  for (int i = 0; i < mNumVtx; i++) {
    mVtxCrd(i,0) = coordinates_element[mNumDim * i + 0];
    mVtxCrd(i,1) = coordinates_element[mNumDim * i + 1];
  }

  // Save element center
  mElmCtr << mVtxCrd.col(0).mean(),
             mVtxCrd.col(1).mean();

}

template <typename Derived>
void QuadNew<Derived>::attachMaterialProperties(const ExodusModel *model, std::string parameter) {
  Vector4d material_at_vertices;
  for (int i = 0; i < mNumVtx; i++) {
    material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(
        mElmCtr, parameter, i);
  }
  mPar[parameter] = material_at_vertices;
}

template <typename Derived>
void QuadNew<Derived>::attachReceiver(std::vector<std::shared_ptr<Receiver>> &receivers) {

  for (auto &rec: receivers) {
    double x1 = rec->PysLocX1();
    double x2 = rec->PysLocX2();
    if (Derived::checkHull(x1, x2, mVtxCrd)) {
      Vector2d ref_loc = Derived::inverseCoordinateTransform(x1, x2, mVtxCrd);
      rec->SetRefLocR(ref_loc(0));
      rec->SetRefLocS(ref_loc(1));
      mRec.push_back(rec);
    }
  }
}

template <typename Derived>
void QuadNew<Derived>::attachSource(std::vector<std::shared_ptr<Source>> sources) {
  for (auto &source: sources) {
    double x1 = source->PhysicalLocationX();
    double x2 = source->PhysicalLocationZ();
    if (Derived::checkHull(x1, x2, mVtxCrd)) {
      Vector2d ref_loc = Derived::inverseCoordinateTransform(x1, x2, mVtxCrd);
      source->setReferenceLocationEps(ref_loc(0));
      source->setReferenceLocationEta(ref_loc(1));
      mSrc.push_back(source);
    }
  }
}

template <typename Derived>
MatrixXd QuadNew<Derived>::setupGradientOperator(const int order) {

  int num_pts_r = QuadNew<Derived>::GllPointsForOrder(order).size();
  int num_pts_s = QuadNew<Derived>::GllPointsForOrder(order).size();
  double eta = QuadNew<Derived>::GllPointsForOrder(order)(0);

  MatrixXd grad(num_pts_s, num_pts_r);
  MatrixXd test(num_pts_s, num_pts_r);
  for (int i = 0; i < num_pts_r; i++) {
    double eps = QuadNew<Derived>::GllPointsForOrder(order)(i);
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
MatrixXd QuadNew<ConcreteShape>::computeGradient(const Ref<const MatrixXd> &field) {

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

      // Optimized gradient for tensorized GLL basis.
      std::tie(invJac, mDetJac(index)) = ConcreteShape::inverseJacobianAtPoint(r, s, mVtxCrd);
      mGradWork.row(index) = invJac * (refGrad <<
        mGrd.row(r_ind).dot(rVectorStride(field, s_ind, mNumIntPtsS, mNumIntPtsR)),
        mGrd.row(s_ind).dot(sVectorStride(field, r_ind, mNumIntPtsS, mNumIntPtsR))).finished();

    }
  }

  return mGradWork;

}

template <typename ConcreteShape>
VectorXd QuadNew<ConcreteShape>::ParAtIntPts(const std::string &par) {

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
VectorXd QuadNew<ConcreteShape>::applyTestAndIntegrate(const Ref<const VectorXd> &f) {

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
VectorXd QuadNew<ConcreteShape>::applyGradTestAndIntegrate(const Ref<const MatrixXd> &f) {

  Matrix2d invJac;
  for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
    for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

      double r = mIntCrdR(r_ind);
      double s = mIntCrdS(s_ind);

      double _;
      std::tie(invJac, _) = ConcreteShape::inverseJacobianAtPoint(r,s,mVtxCrd);

      double dphi_r_dfx = mIntWgtS(s_ind) *
          mIntWgtR.dot(((rVectorStride(mDetJac, s_ind, mNumIntPtsS, mNumIntPtsR)).array() *
          rVectorStride(f.col(0), s_ind, mNumIntPtsS, mNumIntPtsR).array() *
          mGrd.col(r_ind).array()).matrix());

      double dphi_s_dfx = mIntWgtR(r_ind) *
          mIntWgtS.dot(((sVectorStride(mDetJac, r_ind, mNumIntPtsS, mNumIntPtsR)).array() *
          sVectorStride(f.col(0), r_ind, mNumIntPtsS, mNumIntPtsR).array() *
          mGrd.col(s_ind).array()).matrix());

      double dphi_r_dfy = mIntWgtS(s_ind) *
          mIntWgtR.dot(((rVectorStride(mDetJac, s_ind, mNumIntPtsS, mNumIntPtsR)).array() *
          rVectorStride(f.col(1), s_ind, mNumIntPtsS, mNumIntPtsR).array() *
          mGrd.col(r_ind).array()).matrix());

      double dphi_s_dfy = mIntWgtR(r_ind) *
          mIntWgtS.dot(((sVectorStride(mDetJac, r_ind, mNumIntPtsS, mNumIntPtsR)).array() *
          sVectorStride(f.col(1), r_ind, mNumIntPtsS, mNumIntPtsR).array() *
          mGrd.col(s_ind).array()).matrix());

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

template <typename Derived>
double QuadNew<Derived>::integrateField(const Eigen::Ref<const Eigen::VectorXd> &field) {

  double val = 0;
  Matrix2d inverse_Jacobian;
  double detJ;
  for (int i = 0; i < mNumIntPtsS; i++) {
    for (int j = 0; j < mNumIntPtsR; j++) {

      double r = mIntCrdR(j);
      double s = mIntCrdS(i);
      std::tie(inverse_Jacobian, detJ) = Derived::inverseJacobianAtPoint(r, s, mVtxCrd);
      val += field(j + i * mNumIntPtsR) * mIntWgtR(j) *
          mIntWgtS(i) * detJ;

    }
  }

  return val;
}

template <typename Derived>
void QuadNew<Derived>::setBoundaryConditions(Mesh *mesh) {
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

template <typename Derived>
void QuadNew<Derived>::applyDirichletBoundaries(Mesh *mesh, Options &options, const std::string &fieldname) {

  if (! mBndElm) return;

  double value = 0;
  auto dirchlet_boundary_names = options.DirichletBoundaries();
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
template class QuadNew<QuadP1>;

