#include <Mesh/Mesh.h>
#include <Source/Source.h>
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>
#include <Receiver/Receiver.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>

extern "C" {
#include <Element/HyperCube/Autogen/quad_autogen.h>
}

using namespace Eigen;

template<typename ConcreteShape>
TensorQuad<ConcreteShape>::TensorQuad(std::unique_ptr<Options> const &options) {

  /* Ensure we've set parameters correctly. */
  if (options->PolynomialOrder() <= 0 || options->PolynomialOrder() > mMaxOrder) {
    throw std::runtime_error("Polynomial order " + std::to_string(options->PolynomialOrder()) +
    " not supported for quad. Enter a value between 1 and " + std::to_string(mMaxOrder));
  }


  mPlyOrd = options->PolynomialOrder();
  mNumDofVtx = 1;
  mNumDofEdg = mPlyOrd - 1;
  mNumDofFac = (mPlyOrd - 1) * (mPlyOrd - 1);
  mNumDofVol = 0;

  mGrd = TensorQuad<ConcreteShape>::setupGradientOperator(mPlyOrd);
  mIntCrdR = TensorQuad<ConcreteShape>::GllPointsForOrder(mPlyOrd);
  mIntCrdS = TensorQuad<ConcreteShape>::GllPointsForOrder(mPlyOrd);
  mIntWgtR = TensorQuad<ConcreteShape>::GllIntegrationWeightsForOrder(mPlyOrd);
  mIntWgtS = TensorQuad<ConcreteShape>::GllIntegrationWeightsForOrder(mPlyOrd);

  mNumIntPtsS = mIntCrdS.size();
  mNumIntPtsR = mIntWgtR.size();
  mNumIntPnt = mNumIntPtsS * mNumIntPtsR;
  /* Identity closure for tensor basis. */
  mClsMap = IntVec::LinSpaced(mNumIntPnt, 0, mNumIntPnt - 1);

  mDetJac.setZero(mNumIntPnt);
  mParWork.setZero(mNumIntPnt);
  mStiffWork.setZero(mNumIntPnt);
  mGradWork.setZero(mNumIntPnt, mNumDim);

}

template<typename ConcreteShape>
RealVec TensorQuad<ConcreteShape>::GllPointsForOrder(const PetscInt order) {
  if (order > mMaxOrder) { throw std::runtime_error("Polynomial order not supported"); }
  RealVec gll_points(order + 1);
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

template<typename ConcreteShape>
RealVec TensorQuad<ConcreteShape>::GllIntegrationWeightsForOrder(const PetscInt order) {
  if (order > mMaxOrder) { throw std::runtime_error("Polynomial order not supported"); }
  RealVec integration_weights(order + 1);
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

template<typename ConcreteShape>
RealMat TensorQuad<ConcreteShape>::setupGradientOperator(const PetscInt order) {

  if (order > mMaxOrder) { throw std::runtime_error("Polynomial order not supported"); }
  PetscInt num_pts_r = TensorQuad<ConcreteShape>::GllPointsForOrder(order).size();
  PetscInt num_pts_s = TensorQuad<ConcreteShape>::GllPointsForOrder(order).size();
  PetscReal eta = TensorQuad<ConcreteShape>::GllPointsForOrder(order)(0);

  RealMat grad(num_pts_s, num_pts_r);
  RealMat test(num_pts_s, num_pts_r);
  for (PetscInt i = 0; i < num_pts_r; i++) {
    PetscReal eps = TensorQuad<ConcreteShape>::GllPointsForOrder(order)(i);
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

template<typename ConcreteShape>
RealVec TensorQuad<ConcreteShape>::interpolateLagrangePolynomials(const PetscReal r,
                                                                  const PetscReal s,
                                                                  const PetscInt order) {

  if (order > mMaxOrder) { throw std::runtime_error("Polynomial order not supported"); }
  PetscInt n_points = (order + 1) * (order + 1);
  RealVec gll_coeffs(n_points);
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

template<typename ConcreteShape>
RealMat TensorQuad<ConcreteShape>::computeGradient(const Ref<const RealVec> &field) {

  RealVec2 refGrad;
  RealMat2x2 invJac;

  // Loop over all GLL points.
  for (PetscInt s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
    for (PetscInt r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

      // gll index.
      PetscInt index = r_ind + s_ind * mNumIntPtsR;

      // (r,s) coordinates for this point.
      PetscReal r = mIntCrdR(r_ind);
      PetscReal s = mIntCrdS(s_ind);

      // inverse jacobian at this point.
      ConcreteShape::inverseJacobianAtPoint(r, s, mVtxCrd, mDetJac(index), invJac);

      // compute gradient in the reference quad.
      refGrad.setZero(2);
      for (int i = 0; i < mNumIntPtsR; i++) {
        refGrad(0) += mGrd(r_ind, i) * field(i + s_ind * mNumIntPtsR);
        refGrad(1) += mGrd(s_ind, i) * field(r_ind + i * mNumIntPtsR);
      }

      // transform gradient to physical coordinates.
      mGradWork.row(index).noalias() = invJac * refGrad;
    }
  }

  return mGradWork;

}

template<typename ConcreteShape>
void TensorQuad<ConcreteShape>::attachVertexCoordinates(std::unique_ptr<Mesh> const &mesh) {

  Vec coordinates_local;
  PetscInt coordinate_buffer_size;
  PetscSection coordinate_section;
  PetscReal *coordinates_buffer = NULL;

  DMGetCoordinatesLocal(mesh->DistributedMesh(), &coordinates_local);
  DMGetCoordinateSection(mesh->DistributedMesh(), &coordinate_section);
  DMPlexVecGetClosure(mesh->DistributedMesh(), coordinate_section, coordinates_local,
                      mElmNum, &coordinate_buffer_size, &coordinates_buffer);
  std::vector<PetscReal> coordinates_element(coordinates_buffer,
                                             coordinates_buffer + coordinate_buffer_size);
  DMPlexVecRestoreClosure(mesh->DistributedMesh(), coordinate_section, coordinates_local,
                          mElmNum, &coordinate_buffer_size, &coordinates_buffer);

  // Get all coordinates.
  for (PetscInt i = 0; i < mNumVtx; i++) {
    mVtxCrd(i, 0) = coordinates_element[mNumDim * i + 0];
    mVtxCrd(i, 1) = coordinates_element[mNumDim * i + 1];
  }

  // Save edge maps.
  mEdgMap = mesh->EdgeNumbers(mElmNum);

  // Save element center
  mElmCtr << mVtxCrd.col(0).mean(), mVtxCrd.col(1).mean();

}

template<typename ConcreteShape>
void TensorQuad<ConcreteShape>::attachMaterialProperties(std::unique_ptr<ExodusModel> const &model,
                                                         std::string parameter) {
  RealVec4 material_at_vertices;
  for (int i = 0; i < mNumVtx; i++) {
    material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(mElmCtr, parameter, i);
  }
  mPar[parameter] = material_at_vertices;
}

template<typename ConcreteShape>
bool TensorQuad<ConcreteShape>::attachReceiver(std::unique_ptr<Receiver> &receiver,
                                               const bool finalize) {
  if (!receiver) { return false; }
  double x1 = receiver->LocX();
  double x2 = receiver->LocY();
  if (ConcreteShape::checkHull(x1, x2, mVtxCrd)) {
    if (!finalize) { return true; }
    RealVec2 ref_loc = ConcreteShape::inverseCoordinateTransform(x1, x2, mVtxCrd);
    receiver->SetRefLocR(ref_loc(0));
    receiver->SetRefLocS(ref_loc(1));
    mRec.push_back(std::move(receiver));
    return true;
  }
  return false;
}

template<typename ConcreteShape>
bool TensorQuad<ConcreteShape>::attachSource(std::unique_ptr<Source> &source, const bool finalize) {
  if (!source) { return false; }
  PetscReal x1 = source->LocX();
  PetscReal x2 = source->LocY();
  if (ConcreteShape::checkHull(x1, x2, mVtxCrd)) {
    if (!finalize) { return true; }
    RealVec2 ref_loc = ConcreteShape::inverseCoordinateTransform(x1, x2, mVtxCrd);
    source->SetLocR(ref_loc(0));
    source->SetLocS(ref_loc(1));
    mSrc.push_back(std::move(source));
    return true;
  }
  return false;
}

template<typename ConcreteShape>
RealVec TensorQuad<ConcreteShape>::getDeltaFunctionCoefficients(const Eigen::Ref<RealVec>& pnt) {

  PetscReal r = pnt(0), s = pnt(1);
  RealMat2x2 _;
  mParWork = interpolateLagrangePolynomials(r, s, mPlyOrd);
  for (PetscInt s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
    for (PetscInt r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

      PetscReal ri = mIntCrdR(r_ind);
      PetscReal si = mIntCrdS(s_ind);

      PetscReal detJac;
      ConcreteShape::inverseJacobianAtPoint(ri, si, mVtxCrd, detJac, _);

      mParWork(r_ind + s_ind * mNumIntPtsR) /= (mIntWgtR(r_ind) * mIntWgtS(s_ind) * detJac);

    }
  }
  return mParWork;
}

template<typename ConcreteShape>
RealVec TensorQuad<ConcreteShape>::ParAtIntPts(const std::string &par) {

  for (PetscInt s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
    for (PetscInt r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

      PetscReal r = mIntCrdR(r_ind);
      PetscReal s = mIntCrdS(s_ind);
      mParWork(r_ind + s_ind * mNumIntPtsR) =
          ConcreteShape::interpolateAtPoint(r, s).dot(mPar[par]);

    }
  }

  return mParWork;
}

template<typename ConcreteShape>
RealVec TensorQuad<ConcreteShape>::applyTestAndIntegrate(const Ref<const RealVec> &f) {

  RealMat2x2 invJac;
  for (PetscInt s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
    for (PetscInt r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

      // gll index.
      PetscInt index = r_ind + s_ind * mNumIntPtsR;

      // (r,s) coordinate at this point.
      PetscReal r = mIntCrdR(r_ind);
      PetscReal s = mIntCrdS(s_ind);

      PetscReal detJac;
      ConcreteShape::inverseJacobianAtPoint(r, s, mVtxCrd, detJac, invJac);
      mParWork(index) = f(index) * detJac * mIntWgtR(r_ind) * mIntWgtS(s_ind);

    }
  }

  return mParWork;

}

template<typename ConcreteShape>
RealVec TensorQuad<ConcreteShape>::applyGradTestAndIntegrate(const Ref<const RealMat> &f) {

  RealMat2x2 invJac;
  RealVec2 dphi_rs_dfx, dphi_rs_dfy;

  // Loop over all GLL points.
  for (PetscInt s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
    for (PetscInt r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

      // Reset derivatives for this point.
      dphi_rs_dfx.setZero(); dphi_rs_dfy.setZero();

      // Get reference coordinates.
      PetscReal r = mIntCrdR(r_ind);
      PetscReal s = mIntCrdS(s_ind);

      // Loop over the tensor basis. Note we already have detJac at the relevant points.
      for (PetscInt i = 0; i < mNumIntPtsR; i++) {
        PetscInt r_index = i + s_ind * mNumIntPtsR;
        PetscInt s_index = r_ind + i * mNumIntPtsR;
        dphi_rs_dfx(0) += mDetJac(r_index) * f(r_index, 0) * mGrd(i, r_ind) * mIntWgtR(i);
        dphi_rs_dfx(1) += mDetJac(s_index) * f(s_index, 0) * mGrd(i, s_ind) * mIntWgtS(i);
        dphi_rs_dfy(0) += mDetJac(r_index) * f(r_index, 1) * mGrd(i, r_ind) * mIntWgtR(i);
        dphi_rs_dfy(1) += mDetJac(s_index) * f(s_index, 1) * mGrd(i, s_ind) * mIntWgtS(i);
      }

      // Multiply by relevant integration weights.
      dphi_rs_dfx(0) *= mIntWgtS(s_ind);
      dphi_rs_dfx(1) *= mIntWgtR(r_ind);
      dphi_rs_dfy(0) *= mIntWgtS(s_ind);
      dphi_rs_dfy(1) *= mIntWgtR(r_ind);

      // Get the inverse Jacobain again at this point.
      PetscReal _;
      ConcreteShape::inverseJacobianAtPoint(r, s, mVtxCrd, _, invJac);

      // Optimization opportunity.
      // Transform the quantities to physical coordinates.
      PetscInt index = r_ind + s_ind * mNumIntPtsR;
      mStiffWork(index)  = invJac.row(0).dot(dphi_rs_dfx) + invJac.row(1).dot(dphi_rs_dfy);

    }
  }

  return mStiffWork;

}

template<typename ConcreteShape>
RealVec2 TensorQuad<ConcreteShape>::getEdgeNormal(const PetscInt edg) {

  RealVec2 n;
  PetscReal x0, x1, y0, y1;
  for (int i = 0; i < mNumVtx; i++) {
    if (mEdgMap[i] == edg) {
      x0 = mVtxCrd(i, 0);
      x1 = mVtxCrd((i + 1) % mNumVtx, 0);
      y0 = mVtxCrd(i, 1);
      y1 = mVtxCrd((i + 1) % mNumVtx, 1);
    }
  }

  // Should be outwards normal with right hand rule.
  n(0) = +1 * (y1 - y0);
  n(1) = -1 * (x1 - x0);
  return n / n.norm();

}

template <typename ConcreteShape>
void TensorQuad<ConcreteShape>::setEdgeToValue(
    const PetscInt edg, const PetscScalar val, Eigen::Ref<RealVec> f) {

  /* Get proper dofs in the tensor basis. */
  PetscInt start, stride;
  if      (edg == 0) { start = 0;                               stride = 1; }
  else if (edg == 1) { start = mNumIntPtsR - 1;                 stride = mNumIntPtsR; }
  else if (edg == 2) { start = mNumIntPtsR * mNumIntPtsS - 1;   stride = -1; }
  else if (edg == 3) { start = mNumIntPtsR * (mNumIntPtsS - 1); stride = -1 * mNumIntPtsR; }

  for (PetscInt i = start, j = 0; j < mNumIntPtsR; i += stride, j++) {
    f(i) = val;
  }

}


template<typename ConcreteShape>
RealVec TensorQuad<ConcreteShape>::applyTestAndIntegrateEdge(const Eigen::Ref<const RealVec> &f,
                                                             const PetscInt edg) {

  mParWork.setZero();

  // get edge vertices.
  PetscInt start, stride;
  PetscReal x0, x1, y0, y1;
  for (PetscInt i = 0; i < mNumVtx; i++) {
    if (mEdgMap[i] == edg) {
      x0 = mVtxCrd(i, 0);
      x1 = mVtxCrd((i + 1) % mNumVtx, 0);
      y0 = mVtxCrd(i, 1);
      y1 = mVtxCrd((i + 1) % mNumVtx, 1);
      if (i == 0) {
        start = 0;
        stride = 1;
      } else if (i == 1) {
        start = mNumIntPtsR - 1;
        stride = mNumIntPtsR;
      } else if (i == 2) {
        start = mNumIntPtsR * mNumIntPtsS - 1;
        stride = -1;
      } else {
        start = mNumIntPtsR * (mNumIntPtsS - 1);
        stride = -1 * mNumIntPtsR;
      }
    }
  }

  // compute 'edge jacobian'.
  PetscReal dx = x1 - x0;
  PetscReal dy = y1 - y0;
  PetscReal d = sqrt(dx * dx + dy * dy) / 2.0;

  // compute coefficients.
  for (PetscInt i = start, j = 0; j < mNumIntPtsR; i += stride, j++) {
    mParWork(i) = f(i) * d * mIntWgtR(j);
  }

  return mParWork;

}

// Instantiate combinatorical cases.
template
class TensorQuad<QuadP1>;

