#include <Mesh/Mesh.h>
#include <Source/Source.h>
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>
#include <Receiver/Receiver.h>
#include <Element/HyperCube/HexP1.h>
#include <Element/HyperCube/Hexahedra.h>

#include <complex>

// Extern.
extern "C" {
#include <Element/HyperCube/Autogen/quad_autogen.h>
#include <Element/HyperCube/Autogen/hex_autogen.h>
}

using namespace Eigen;

template <typename ConcreteHex>
Hexahedra<ConcreteHex>::Hexahedra(std::unique_ptr<Options> const &options) {

  /* Ensure we've set parameters correctly. */
  assert(options->PolynomialOrder() > 0);

  // Basic properties.
  mPlyOrd = options->PolynomialOrder();
  
  // Gll points.
  mNumDofVtx = 1;
  mNumDofEdg = mPlyOrd - 1;
  mNumDofFac = (mPlyOrd - 1) * (mPlyOrd - 1);
  mNumDofVol = (mPlyOrd - 1) * (mPlyOrd - 1) * (mPlyOrd - 1);

  // Integration points.
  mIntCrdR = GllPointsForOrder(mPlyOrd);
  mIntCrdS = GllPointsForOrder(mPlyOrd);
  mIntCrdT = GllPointsForOrder(mPlyOrd);
  mIntWgtR = GllIntegrationWeights(mPlyOrd);
  mIntWgtS = GllIntegrationWeights(mPlyOrd);
  mIntWgtT = GllIntegrationWeights(mPlyOrd);
  
  // Save number of integration points.
  mNumIntPtsR = mIntCrdR.size();
  mNumIntPtsS = mIntCrdS.size();
  mNumIntPtsT = mIntCrdT.size();
  mNumIntPnt = mNumIntPtsR * mNumIntPtsS * mNumIntPtsT;

  // setup evaluated derivatives of test functions
  mGrd = setupGradientOperator(mPlyOrd);
  mGrdT = mGrd.transpose();

  mGrdWgt.resize(mNumIntPtsR,mNumIntPtsR);
  for(PetscInt i=0;i<mNumIntPtsR;i++) {
    for(PetscInt j=0;j<mNumIntPtsR;j++) {
      mGrdWgt(i,j) = mGrd(i,j)*mIntWgtR[i];
    }
  }
  mGrdWgtT = mGrdWgt.transpose();
  
  mDetJac.setZero(mNumIntPnt);
  mParWork.setZero(mNumIntPnt);
  mStiffWork.setZero(mNumIntPnt);
  mGradWork.setZero(mNumIntPnt, mNumDim);

}

template <typename ConcreteHex>
RealVec Hexahedra<ConcreteHex>::GllPointsForOrder(const PetscInt order) {
  RealVec gll_points(order + 1);
  if (order == 1) {
    gll_coordinates_order1_square(gll_points.data());
  } else if (order == 2) {
    gll_coordinates_order2_square(gll_points.data());
  } else if (order == 3) {
    gll_coordinates_order3_square(gll_points.data());
  }
  return gll_points;
}

template <typename ConcreteHex>
RealVec Hexahedra<ConcreteHex>::GllIntegrationWeights(const PetscInt order) {
  RealVec integration_weights(order + 1);
  if (order == 1) {
    gll_weights_order1_square(integration_weights.data());
  } else if (order == 2) {
    gll_weights_order2_square(integration_weights.data());
  } else if (order == 3) {
    gll_weights_order3_square(integration_weights.data());
  }
  return integration_weights;
}

template <typename ConcreteHex>
void Hexahedra<ConcreteHex>::attachVertexCoordinates(std::unique_ptr<Mesh> const &mesh) {

  Vec coordinates_local;
  PetscInt coordinate_buffer_size;
  PetscSection coordinate_section;
  PetscReal *coordinates_buffer = NULL;


  DMGetCoordinatesLocal(mesh->DistributedMesh(), &coordinates_local);
  DMGetCoordinateSection(mesh->DistributedMesh(), &coordinate_section);
  DMPlexVecGetClosure(mesh->DistributedMesh(), coordinate_section, coordinates_local, mElmNum,
                      &coordinate_buffer_size, &coordinates_buffer);
  std::vector<PetscReal> coordinates_element(coordinates_buffer,
                                             coordinates_buffer + coordinate_buffer_size);
  DMPlexVecRestoreClosure(mesh->DistributedMesh(), coordinate_section, coordinates_local, mElmNum,
                          &coordinate_buffer_size, &coordinates_buffer);

  for (PetscInt i = 0; i < mNumVtx; i++) {
    mVtxCrd(i, 0) = coordinates_element[mNumDim * i + 0];
    mVtxCrd(i, 1) = coordinates_element[mNumDim * i + 1];
    mVtxCrd(i, 2) = coordinates_element[mNumDim * i + 2];
  }

  // Save element center
  mElmCtr <<
    mVtxCrd.col(0).mean(),
    mVtxCrd.col(1).mean(),
    mVtxCrd.col(2).mean();

}

template <typename ConcreteHex>
void Hexahedra<ConcreteHex>::setEdgeToValue(
    const PetscInt edg, const PetscReal val, Eigen::Ref<RealVec> f) {

  /* Get proper dofs in the tensor basis. */
  PetscInt s_sta, s_str, r_sta, r_str;
  switch (edg) {
    case 0: /* bottom */
      s_sta = 0;
      r_sta = 0;
      s_str = 1;
      r_str = 1;
      break;
    case 1: /* top */
      s_sta = 0;
      r_sta = mNumIntPnt - mNumIntPtsR * mNumIntPtsS;
      s_str = 1;
      r_str = 1;
      break;
    case 2: /* left */
      s_sta = 0;
      r_sta = 0;
      r_str = 1;
      s_str = mNumIntPtsS;
      break;
    case 3: /* right */
      r_sta = mNumIntPtsR * mNumIntPtsS - mNumIntPtsS;
      s_sta = 0;
      r_str = 1;
      s_str = mNumIntPtsR;
      break;
    case 4: /* back */
      r_sta = mNumIntPtsR - 1;
      s_sta = 0;
      r_str = mNumIntPtsS;
      s_str = mNumIntPtsR;
      break;
    case 5: /* front */
      r_sta = 0;
      s_sta = 0;
      r_str = mNumIntPtsS;
      s_str = mNumIntPtsR;
      break;
  }

  /* Set the proper dofs. */
  for (PetscInt s = s_sta, icnt = 0; icnt < mNumIntPtsR; s += s_str, icnt++) {
    for (PetscInt r = r_sta, jcnt = 0; jcnt < mNumIntPtsS; r += r_str, jcnt++) {
      f(r + s * mNumIntPtsS) = val;
    }
  }

}

template <typename ConcreteHex>
RealVec Hexahedra<ConcreteHex>::applyTestAndIntegrateEdge(const Eigen::Ref<const RealVec> &f,
                                                          const PetscInt edg) {

  PetscInt s_sta, r_sta, s_str, r_str;
  RealVec3 q0, q1, q2, q3;

  /* Get vertices defining edge. */
  switch (edg) {
    case 0: /* bottom */
      q0 = mVtxCrd.row(0); q1 = mVtxCrd.row(1);
      q2 = mVtxCrd.row(2); q3 = mVtxCrd.row(3);
      s_sta = 0;
      r_sta = 0;
      s_str = 1;
      r_str = 1;
      break;
    case 1: /* top */
      q0 = mVtxCrd.row(4); q1 = mVtxCrd.row(5);
      q2 = mVtxCrd.row(6); q3 = mVtxCrd.row(7);
      s_sta = 0;
      r_sta = mNumIntPnt - mNumIntPtsR * mNumIntPtsS;
      s_str = 1;
      r_str = 1;
      break;
    case 2: /* left */
      q0 = mVtxCrd.row(0); q1 = mVtxCrd.row(4);
      q2 = mVtxCrd.row(7); q3 = mVtxCrd.row(1);
      s_sta = 0;
      r_sta = 0;
      r_str = 1;
      s_str = mNumIntPtsS;
      break;
    case 3: /* right */
      q0 = mVtxCrd.row(3); q1 = mVtxCrd.row(2);
      q2 = mVtxCrd.row(6); q3 = mVtxCrd.row(5);
      r_sta = mNumIntPtsR * mNumIntPtsS - mNumIntPtsS;
      s_sta = 0;
      r_str = 1;
      s_str = mNumIntPtsR;
      break;
    case 4: /* back */
      q0 = mVtxCrd.row(1); q1 = mVtxCrd.row(7);
      q2 = mVtxCrd.row(6); q3 = mVtxCrd.row(2);
      r_sta = mNumIntPtsR - 1;
      s_sta = 0;
      r_str = mNumIntPtsS;
      s_str = mNumIntPtsR;
      break;
    case 5: /* front */
      q0 = mVtxCrd.row(0); q1 = mVtxCrd.row(3);
      q2 = mVtxCrd.row(5); q3 = mVtxCrd.row(4);
      r_sta = 0;
      s_sta = 0;
      r_str = mNumIntPtsS;
      s_str = mNumIntPtsR;
      break;
  }

  /* Vectors defining face plane. */
  RealVec3 v0 = q1 - q0;
  RealVec3 v1 = q3 - q0;

  /* Face normal. */
  RealVec3 n = v0.cross(v1).normalized();

  /* Compute orthogonal 2D coordinates. */
  RealVec3 e0 = v0.normalized();
  RealVec3 e1 = n.cross(v0).normalized();

  /* Compute coordinates in 2D planar coordinates. */
  QuadVtx eVtx;
  eVtx(0, 0) = e0.dot(q0 - q0);
  eVtx(1, 0) = e0.dot(q1 - q0);
  eVtx(2, 0) = e0.dot(q2 - q0);
  eVtx(3, 0) = e0.dot(q3 - q0);
  eVtx(0, 1) = e1.dot(q0 - q0);
  eVtx(1, 1) = e1.dot(q1 - q0);
  eVtx(2, 1) = e1.dot(q2 - q0);
  eVtx(3, 1) = e1.dot(q3 - q0);

  mParWork.setZero();
  for (PetscInt s_ind = s_sta, s_cnt = 0; s_cnt < mNumIntPtsR; s_ind += s_str, s_cnt++) {
    for (PetscInt r_ind = r_sta, r_cnt = 0; r_cnt < mNumIntPtsS; r_ind += r_str, r_cnt++) {

      PetscReal detJac;
      PetscReal r = mIntCrdR(r_cnt);
      PetscReal s = mIntCrdS(s_cnt);
      PetscInt ind = r_ind + s_ind * mNumIntPtsS;

      ConcreteHex::faceJacobianAtPoint(r, s, eVtx, detJac);
      mParWork(ind) = f(ind) * detJac * mIntWgtR(r_cnt) * mIntWgtS(s_cnt);

    }
  }

  return mParWork;

}


template <typename ConcreteHex>
void Hexahedra<ConcreteHex>::attachMaterialProperties(
    std::unique_ptr<ExodusModel> const &model, std::string parameter_name) {

  RealVec material_at_vertices(mNumVtx);

  for (auto i = 0; i < mNumVtx; i++) {
    material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(mElmCtr, parameter_name, i);
  }
  mPar[parameter_name] = material_at_vertices;
  
}



template <typename ConcreteHex>
bool Hexahedra<ConcreteHex>::attachReceiver(std::unique_ptr<Receiver> &receiver,
                                            const bool finalize) {
  if (!receiver) { return false; }
  PetscReal x1 = receiver->LocX();
  PetscReal x2 = receiver->LocY();
  PetscReal x3 = receiver->LocZ();
  if (ConcreteHex::checkHull(x1, x2, x3, mVtxCrd)) {
    if (!finalize) { return true; }
    RealVec3 ref_loc = ConcreteHex::inverseCoordinateTransform(x1, x2, x3, mVtxCrd);
    receiver->SetRefLocR(ref_loc(0));
    receiver->SetRefLocS(ref_loc(1));
    receiver->SetRefLocT(ref_loc(2));
    mRec.push_back(std::move(receiver));
    return true;
  }
  return false;
}

template <typename ConcreteHex>
bool Hexahedra<ConcreteHex>::attachSource(std::unique_ptr<Source> &source, const bool finalize) {
  if (!source) { return false; }
  PetscReal x1 = source->LocX();
  PetscReal x2 = source->LocY();
  PetscReal x3 = source->LocZ();
  if (ConcreteHex::checkHull(x1, x2, x3, mVtxCrd)) {
    if (!finalize) { return true; }
    RealVec3 ref_loc = ConcreteHex::inverseCoordinateTransform(x1, x2, x3, mVtxCrd);
    source->SetLocR(ref_loc(0));
    source->SetLocS(ref_loc(1));
    source->SetLocT(ref_loc(2));
    mSrc.push_back(std::move(source));
    return true;
  }
  return false;
}

template <typename ConcreteHex>
RealVec Hexahedra<ConcreteHex>::getDeltaFunctionCoefficients(const Eigen::Ref<RealVec>& pnt) {

  PetscReal r = pnt(0), s = pnt(1), t = pnt(2);
  RealMat3x3 invJ;
  mParWork = interpolateLagrangePolynomials(r, s, t, mPlyOrd);
  for (PetscInt t_ind = 0; t_ind < mNumIntPtsT; t_ind++) {
    for (PetscInt s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (PetscInt r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {


        PetscReal ri = mIntCrdR(r_ind);
        PetscReal si = mIntCrdS(s_ind);
        PetscReal ti = mIntCrdT(t_ind);

        PetscReal detJac;
        ConcreteHex::inverseJacobianAtPoint(ri, si, ti, mVtxCrd, detJac, invJ);
        mParWork(r_ind + s_ind * mNumIntPtsR + t_ind*mNumIntPtsR*mNumIntPtsS) /=
          (mIntWgtR(r_ind) * mIntWgtS(s_ind) * mIntWgtT(t_ind)  * detJac);

      }
    }
  }
  return mParWork;
}

  
template <typename ConcreteHex>
RealVec Hexahedra<ConcreteHex>::interpolateLagrangePolynomials(const PetscReal r,
                                                                   const PetscReal s,
                                                                   const PetscReal t,
                                                                   const PetscInt order)

{

  PetscInt n_points = (order + 1) * (order + 1) * (order + 1);
  RealVec gll_coeffs(n_points);
  if (order == 1) {
    interpolate_order1_hex(r, s, t, gll_coeffs.data());
  } else if (order == 2) {
    interpolate_order2_hex(r, s, t, gll_coeffs.data());
  } else if (order == 3) {
    interpolate_order3_hex(r, s, t, gll_coeffs.data());
  }

  return gll_coeffs;
  
}

template <typename ConcreteHex>
RealMat  Hexahedra<ConcreteHex>::setupGradientOperator(const PetscInt order) {

  auto rn = GllPointsForOrder(order);
  PetscInt num_pts_r = rn.size();
  PetscInt num_pts_s = rn.size();
  PetscReal r = rn(0);
  PetscReal s = rn(0);
  
  RealMat  grad(num_pts_s, num_pts_r);
  RealMat test(num_pts_r, num_pts_s);
  for (PetscInt i=0; i<num_pts_r; i++) {
    PetscReal r = rn[i];
    if (order == 1) {
      interpolate_eps_derivative_order1_square(s, test.data());
    } else if (order == 2) {
      interpolate_eps_derivative_order2_square(r, s, test.data());
    } else if (order == 3) {
      interpolate_eps_derivative_order3_square(r, s, test.data());
    }
    grad.row(i) = test.col(0);
  }
  
  return grad;
}

template <typename ConcreteHex>
RealMat  Hexahedra<ConcreteHex>::computeGradient(const Ref<const RealVec> &field) {

  RealMat3x3 invJac;
  RealVec3 refGrad;

  // Loop over all GLL points.
  for (PetscInt t_ind = 0; t_ind < mNumIntPtsT; t_ind++) {
    for (PetscInt s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (PetscInt r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        // gll index.
        PetscInt index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;

        // (r,s,t) coordinates for this point.
        PetscReal r = mIntCrdR(r_ind);
        PetscReal s = mIntCrdS(s_ind);
        PetscReal t = mIntCrdT(t_ind);

        // Optimized gradient for tensorized GLL basis.
        PetscReal detJ;
        ConcreteHex::inverseJacobianAtPoint(r, s, t, mVtxCrd, detJ, invJac);
        
        // mGradWork.row(index) = invJac * (refGrad <<
        //                                  mGrd.row(r_ind).dot(rVectorStride(field,s_ind,t_ind,
        //                                                                    mNumIntPtsR,mNumIntPtsS,mNumIntPtsT)),
        //                                  mGrd.row(s_ind).dot(sVectorStride(field, r_ind, t_ind,
        //                                                                    mNumIntPtsR,mNumIntPtsS,mNumIntPtsT)),
        //                                  mGrd.row(t_ind).dot(tVectorStride(field, r_ind, s_ind,
        //                                                                    mNumIntPtsR,mNumIntPtsS,mNumIntPtsT))).finished();
        
        // loop by hand is 10 us faster (27 -> 17)
        refGrad.setZero(3);
        for(PetscInt i=0;i<mNumIntPtsR;i++) {
          refGrad(0) += mGrd(r_ind,i)*field(i + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS);
          refGrad(1) += mGrd(s_ind,i)*field(r_ind + i * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS);
          refGrad(2) += mGrd(t_ind,i)*field(r_ind + s_ind * mNumIntPtsR + i * mNumIntPtsR * mNumIntPtsS);
        }
        mGradWork.row(index) = invJac * refGrad;

      }
    }
  }

  return mGradWork;

}

template <typename ConcreteHex>
RealVec Hexahedra<ConcreteHex>::ParAtIntPts(const std::string &par) {

  // Loop over all GLL points.
  for (PetscInt t_ind = 0; t_ind < mNumIntPtsT; t_ind++) {
    for (PetscInt s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (PetscInt r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        // gll index.
        PetscInt index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;

        // (r,s,t) coordinates for this point.
        PetscReal r = mIntCrdR(r_ind);
        PetscReal s = mIntCrdS(s_ind);
        PetscReal t = mIntCrdT(t_ind);

        mParWork(index) = ConcreteHex::interpolateAtPoint(r,s,t).dot(mPar[par]);
      }
    }
  }

  return mParWork;
}

template <typename ConcreteHex>
RealVec Hexahedra<ConcreteHex>::applyTestAndIntegrate(const Ref<const RealVec> &f) {

  PetscInt i = 0;
  PetscReal detJac;
  RealMat3x3 invJac;
  RealVec result(mNumIntPnt);
  for (PetscInt t_ind = 0; t_ind < mNumIntPtsS; t_ind++) {
    for (PetscInt s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (PetscInt r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        // gll index.
        PetscInt index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;
        
        // (r,s) coordinate at this point.
        PetscReal r = mIntCrdR(r_ind);
        PetscReal s = mIntCrdS(s_ind);
        PetscReal t = mIntCrdT(t_ind);

        ConcreteHex::inverseJacobianAtPoint(r,s,t,mVtxCrd, detJac, invJac);
        result(index) = f(index) * detJac * mIntWgtR(r_ind) * mIntWgtS(s_ind) * mIntWgtS(t_ind);

      }
    }
  }

  return result;

}

template <typename ConcreteHex>
RealVec Hexahedra<ConcreteHex>::applyGradTestAndIntegrate(const Ref<const RealMat > &f) {

  // computes the rotatation into x-y-z, which would normally happen later with more terms.
  RealMat3x3 invJac;
  RealVec3 fi;
  RealMat  fxyz(f.rows(),3);
  for (PetscInt t_ind = 0; t_ind < mNumIntPtsT; t_ind++) {
    for (PetscInt s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (PetscInt r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        // gll index.
        PetscReal r = mIntCrdR(r_ind), s = mIntCrdS(s_ind), t = mIntCrdT(t_ind);
        PetscInt index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;
        ConcreteHex::inverseJacobianAtPoint(r, s, t, mVtxCrd, mDetJac(index), invJac);
        fi << f(index,0),f(index,1),f(index,2);
        fi = invJac.transpose()*fi;

        fxyz(index,0) = fi[0];
        fxyz(index,1) = fi[1];
        fxyz(index,2) = fi[2];
        
      }
    }
  }


  for (PetscInt t_ind = 0; t_ind < mNumIntPtsT; t_ind++) {
    for (PetscInt s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (PetscInt r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        // map reference gradient (lr,ls,lt) to this element (lx,ly,lz)
        auto lr = mGrd.col(r_ind);
        auto ls = mGrd.col(s_ind);
        auto lt = mGrd.col(t_ind);

        PetscInt index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;

        PetscReal dphi_r_dfx = 0;
        PetscReal dphi_s_dfy = 0;
        PetscReal dphi_t_dfz = 0;

        for(PetscInt i=0;i<mNumIntPtsR;i++) {
          
          PetscInt r_index = i + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;
          PetscInt s_index = r_ind + i * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;
          PetscInt t_index = r_ind + s_ind * mNumIntPtsR + i * mNumIntPtsR * mNumIntPtsS;
          
          dphi_r_dfx += mDetJac[r_index] * fxyz(r_index,0) * lr[i] * mIntWgtR[i];
          dphi_s_dfy += mDetJac[s_index] * fxyz(s_index,1) * ls[i] * mIntWgtR[i];
          dphi_t_dfz += mDetJac[t_index] * fxyz(t_index,2) * lt[i] * mIntWgtR[i];
        
        }
        dphi_r_dfx *= mIntWgtR(s_ind) * mIntWgtR(t_ind);
        dphi_s_dfy *= mIntWgtR(r_ind) * mIntWgtR(t_ind);
        dphi_t_dfz *= mIntWgtR(r_ind) * mIntWgtR(s_ind);
        
        mStiffWork(index) = dphi_r_dfx + dphi_s_dfy + dphi_t_dfz;
        
                
      }
    }
  }

  return mStiffWork;

}

//template <typename ConcreteHex>
//void Hexahedra<ConcreteHex>::precomputeConstants() {
//
//  PetscInt num_pts = mNumIntPtsR*mNumIntPtsS*mNumIntPtsT;
//  mDetJac.resize(num_pts);
//  mInvJac.resize(num_pts);
//  RealMat3x3 invJac;
//  // Loop over all GLL points.
//  for (PetscInt t_ind = 0; t_ind < mNumIntPtsT; t_ind++) {
//    for (PetscInt s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
//      for (PetscInt r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {
//
//        // gll index.
//        PetscInt index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;
//
//        // (r,s,t) coordinates for this point.
//        PetscReal r = mIntCrdR(r_ind);
//        PetscReal s = mIntCrdS(s_ind);
//        PetscReal t = mIntCrdT(t_ind);
//
//        // Optimized gradient for tensorized GLL basis.
//        ConcreteHex::inverseJacobianAtPoint(r, s, t, mVtxCrd, mDetJac(index), invJac);
//        mInvJac[index] = invJac;
//      }
//    }
//  }
//}

// Instantiate base case.
template class Hexahedra<HexP1>;

