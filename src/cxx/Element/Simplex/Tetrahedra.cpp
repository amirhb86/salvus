#include <mpi.h>
#include <Mesh/Mesh.h>
#include <Source/Source.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Utilities/Logging.h>
#include <Element/Simplex/TetP1.h>
#include <Element/Simplex/Tetrahedra.h>

extern "C" {
#include <Element/Simplex/Autogen/p3_tetrahedra.h>
}

using namespace Eigen;

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  out << "[";
  size_t last = v.size() - 1;
  for(size_t i = 0; i < v.size(); ++i) {
    out << v[i];
    if (i != last) 
      out << ", ";
  }
  out << "]";
  return out;
}

template <typename ConcreteShape>
MatrixXd Tetrahedra<ConcreteShape>::mGradientOperator;
template <typename ConcreteShape>
VectorXd Tetrahedra<ConcreteShape>::mIntegrationWeights;
template <typename ConcreteShape>
VectorXd Tetrahedra<ConcreteShape>::mIntegrationCoordinates_r;
template <typename ConcreteShape>
VectorXd Tetrahedra<ConcreteShape>::mIntegrationCoordinates_s;
template <typename ConcreteShape>
VectorXd Tetrahedra<ConcreteShape>::mIntegrationCoordinates_t;
template <typename ConcreteShape>
MatrixXd Tetrahedra<ConcreteShape>::mGradientPhi_dr;
template <typename ConcreteShape>
MatrixXd Tetrahedra<ConcreteShape>::mGradientPhi_ds;
template <typename ConcreteShape>
MatrixXd Tetrahedra<ConcreteShape>::mGradientPhi_dt;
template <typename ConcreteShape>
Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Tetrahedra<ConcreteShape>::mGradientPhi_dr_t;
template <typename ConcreteShape>
Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Tetrahedra<ConcreteShape>::mGradientPhi_ds_t;
template <typename ConcreteShape>
Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Tetrahedra<ConcreteShape>::mGradientPhi_dt_t;

// Gets vertex ids from PETSc element closure
std::vector<int> getVertsFromPoint(int point, int numVerts, DM &distributed_mesh) {
  int* points = NULL;
  int numPoints;  
  DMPlexGetTransitiveClosure(distributed_mesh,point,PETSC_TRUE,&numPoints,&points);
  std::vector<int> verts(numVerts);
  // vertices are the last `numVerts` entries in points
  for(int i=numPoints-numVerts;i<numPoints;i++) {
    verts[i-(numPoints-numVerts)] = points[2*i];
  }
  return verts;
}

template <typename ConcreteShape>
std::tuple<VectorXd,VectorXd,VectorXd> Tetrahedra<ConcreteShape>::QuadraturePoints(const int order) {

  VectorXd rn, sn, tn;
  if(order == 3) {
    int num_pts = 50;
    rn.setZero(num_pts);
    sn.setZero(num_pts);
    tn.setZero(num_pts);
    coordinates_p3_tetrahedra_rn(rn.data());
    coordinates_p3_tetrahedra_sn(sn.data());
    coordinates_p3_tetrahedra_tn(tn.data());
  }
  else {
    std::cerr << "ERROR: Order NOT implemented!...\n";
    MPI::COMM_WORLD.Abort(-1);
  }
  return std::make_tuple(rn,sn,tn);

}

template <typename ConcreteShape>
VectorXd Tetrahedra<ConcreteShape>::QuadratureIntegrationWeights(const int order) {

  VectorXd wn;
  if(order == 3) {
    int num_pts = 50;
    wn.setZero(num_pts);
    quadrature_weights_p3_tetrahedra(wn.data());
  } else {
    std::cerr << "ERROR: Order NOT implemented!\n";
    MPI::COMM_WORLD.Abort(-1);
  }
  return wn;
}

template <typename ConcreteShape>
VectorXi Tetrahedra<ConcreteShape>::ClosureMapping(const int order, const int dimension,
                                                   DM &distributed_mesh) {

  VectorXi linear_mapping;
  if(order == 3) {

    linear_mapping.resize(50);
    for(int i=0;i<50;i++) { linear_mapping[i] = i; }

  } else {
    ERROR() << "Order " << order << " tetrahedra closure mapping not implemented!";
  }
  return linear_mapping;
  
}

template <typename ConcreteShape>
Tetrahedra<ConcreteShape>::Tetrahedra(std::unique_ptr<Options> const &options) {

  // Basic properties.
  mPlyOrd = options->PolynomialOrder();
  if(mPlyOrd == 3) {
    // total number of nodes
    mNumIntPnt = 50;
    // degrees of freedom per entity (vertex,edge,face,interior)
    mNumDofVtx = 1;
    mNumDofEdg = 2;
    mNumDofFac = 6;
    mNumDofVol = 10;
  }
  else {
    std::cerr << "ERROR: Order NOT implemented!\n";
    MPI::COMM_WORLD.Abort(-1);
  }
  
  // Integration points and weights
  std::tie(mIntegrationCoordinates_r,mIntegrationCoordinates_s,mIntegrationCoordinates_t) =
    Tetrahedra<ConcreteShape>::QuadraturePoints(options->PolynomialOrder());
  mIntegrationWeights = Tetrahedra<ConcreteShape>::QuadratureIntegrationWeights(options->PolynomialOrder());
          
  setupGradientOperator();

  mParWork.setZero(mNumIntPnt);
  mStiffWork.setZero(mNumIntPnt);
  mGradWork.setZero(mNumIntPnt, mNumDim);
  mGrad_r.resize(mNumIntPnt);
  mGrad_s.resize(mNumIntPnt);
  mGrad_t.resize(mNumIntPnt);
  
}

template <typename ConcreteShape>
void Tetrahedra<ConcreteShape>::attachVertexCoordinates(std::unique_ptr<Mesh> const &mesh) {

  mClsMap = ClosureMapping(3, 3, mesh->DistributedMesh());

  Vec coordinates_local;
  PetscInt coordinate_buffer_size;
  PetscSection coordinate_section;
  PetscReal *coordinates_buffer = NULL;

  DMGetCoordinatesLocal(mesh->DistributedMesh(), &coordinates_local);
  DMGetCoordinateSection(mesh->DistributedMesh(), &coordinate_section);
  DMPlexVecGetClosure(mesh->DistributedMesh(), coordinate_section, coordinates_local, mElmNum,
                      &coordinate_buffer_size, &coordinates_buffer);
  std::vector<PetscReal> coordinates_element(coordinates_buffer, coordinates_buffer+coordinate_buffer_size);
  DMPlexVecRestoreClosure(mesh->DistributedMesh(), coordinate_section, coordinates_local, mElmNum,
                          &coordinate_buffer_size, &coordinates_buffer);
    
  for (int i = 0; i < mNumVtx; i++) {        
    mVtxCrd(i,0) = coordinates_element[mNumDim*i+0];
    mVtxCrd(i,1) = coordinates_element[mNumDim*i+1];
    mVtxCrd(i,2) = coordinates_element[mNumDim*i+2];
  }

  // Save element center
  mElmCtr <<
    mVtxCrd.col(0).mean(),
    mVtxCrd.col(1).mean(),
    mVtxCrd.col(2).mean();

}

template <typename ConcreteShape>
VectorXd Tetrahedra<ConcreteShape>::ParAtIntPts(const std::string &par) {

  // interpolate velocity at all nodes
  for(int i=0;i<mNumIntPnt;i++) {
    auto r = mIntegrationCoordinates_r[i];
    auto s = mIntegrationCoordinates_s[i];
    auto t = mIntegrationCoordinates_t[i];
    mParWork(i) = ConcreteShape::interpolateAtPoint(r,s,t).dot(mPar[par]);
  }  
  
  return mParWork;
}

template <typename ConcreteShape>
void Tetrahedra<ConcreteShape>::attachMaterialProperties(std::unique_ptr<ExodusModel> const &model, std::string parameter_name) {

  Vector4d material_at_vertices;

  for (auto i = 0; i < mNumVtx; i++) {
    material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(mElmCtr, parameter_name, i);
  }
  
  mPar[parameter_name] = material_at_vertices;

}

template <typename ConcreteShape>
bool Tetrahedra<ConcreteShape>::attachSource(std::unique_ptr<Source> &source, const bool finalize) {

  if (!source) { return false; }
  if (ConcreteShape::checkHull(source->LocX(),
                               source->LocY(),
                               source->LocZ(),
                               mVtxCrd)) {
    if (!finalize) { return true; }
    Vector3d reference_location = ConcreteShape::inverseCoordinateTransform(source->LocX(),
                                                                            source->LocY(),
                                                                            source->LocZ(),
                                                                            mVtxCrd);
    source->SetLocR(reference_location(0));
    source->SetLocS(reference_location(1));
    source->SetLocT(reference_location(2));
    mSrc.push_back(std::move(source));
    return true;
  }
  return false;
}

template <typename ConcreteShape>
MatrixXd Tetrahedra<ConcreteShape>::computeGradient(const Ref<const VectorXd>& field) {

  Vector3d phyGrad;
  Vector3d refGrad;

  mGradWork.col(0) = mGradientPhi_dx*field;
  mGradWork.col(1) = mGradientPhi_dy*field;
  mGradWork.col(2) = mGradientPhi_dz*field;
  
  return mGradWork;
}

template <typename ConcreteShape>
VectorXd Tetrahedra<ConcreteShape>::applyGradTestAndIntegrate(const Ref<const MatrixXd>& f) {
  Vector3d refGrad;
  Vector3d phyGrad;

  for(int i=0;i<mNumIntPnt;i++) {
    refGrad << f(i,0), f(i,1), f(i,2);
    phyGrad = mInvJacT * refGrad;
    mGradWork(i,0) = phyGrad(0);
    mGradWork(i,1) = phyGrad(1);
    mGradWork(i,2) = phyGrad(2);
  }
  
  for(int i=0;i<mNumIntPnt;i++) {

    mStiffWork(i) = 0;
    for(int j=0;j<mNumIntPnt;j++) {
      mStiffWork(i) +=
        mIntegrationWeights[j] * (mGradWork(j,0) * mGradientPhi_dr(i,j) +
                                  mGradWork(j,1) * mGradientPhi_ds(i,j) +
                                  mGradWork(j,2) * mGradientPhi_dt(i,j));
    }

    mStiffWork(i) *= mDetJac;
  }
  return mStiffWork;
  
}

template <typename ConcreteShape>
VectorXd Tetrahedra<ConcreteShape>::computeStiffnessFull(const Ref<const VectorXd>& field,
                                                         const Ref<const VectorXd>& vp2) {

  Vector3d refGrad;
  Vector3d phyGrad;
  int num_pts = field.size();
  
  // mGrad_r = mGradientPhi_dr_t * field;
  // mGrad_s = mGradientPhi_ds_t * field;
  // mGrad_t = mGradientPhi_dt_t * field;
  
  for(int i=0;i<mNumIntPnt;i++) {

    // eigen
    // refGrad(0) = mGradientPhi_dr.col(i).dot(field);
    // refGrad(1) = mGradientPhi_ds.col(i).dot(field);
    // refGrad(2) = mGradientPhi_dt.col(i).dot(field);
    
    // no eigen
    refGrad.setZero(3);
    for(int j=0;j<mNumIntPnt;j++) {
      refGrad(0) += mGradientPhi_dr(j,i) * field(j);
      refGrad(1) += mGradientPhi_ds(j,i) * field(j);
      refGrad(2) += mGradientPhi_dt(j,i) * field(j);
    }
    
    // refGrad << mGrad_r(i), mGrad_s(i), mGrad_t(i);
    
    phyGrad = (mInvJacT_x_invJac) * refGrad;
    
    mGradWork(i,0) = vp2(i)*phyGrad(0);
    mGradWork(i,1) = vp2(i)*phyGrad(1);
    mGradWork(i,2) = vp2(i)*phyGrad(2);
    
  }

  for(int i=0;i<mNumIntPnt;i++) {

    refGrad.setZero(3);
    mStiffWork(i) = 0;
    for(int j=0;j<mNumIntPnt;j++) {
      mStiffWork(i) +=
        mIntegrationWeights[j] * (mGradWork(j,0) * mGradientPhi_dr(i,j) +
                                  mGradWork(j,1) * mGradientPhi_ds(i,j) +
                                  mGradWork(j,2) * mGradientPhi_dt(i,j));
    }

    mStiffWork(i) *= mDetJac;
    
    // // eigen version (10-20% slower)
    // mStiffWork(i) = mDetJac * mIntegrationWeights.dot( (fGrad.col(0).array() * mGradientPhi_dr_t.col(i).array() +
    //                                        fGrad.col(1).array() * mGradientPhi_ds_t.col(i).array() +
    //                                        fGrad.col(2).array() * mGradientPhi_dt_t.col(i).array()).matrix());

    
  }
  return mStiffWork;
}

template <typename ConcreteShape>
void Tetrahedra<ConcreteShape>::precomputeElementTerms()
 {
   std::tie(mInvJac,mDetJac) = ConcreteShape::inverseJacobian(mVtxCrd);
   mInvJacT = mInvJac.transpose();
   mElementStiffnessMatrix = buildStiffnessMatrix(ParAtIntPts("VP"));    
 }

template <typename ConcreteShape>
MatrixXd Tetrahedra<ConcreteShape>::buildStiffnessMatrix(const Ref<const VectorXd>& vp2) {
  
  Eigen::Matrix3d invJ;
  double detJ;
  std::tie(invJ,detJ) = ConcreteShape::inverseJacobian(mVtxCrd);

  // mGradientPhi_dr_t = mGradientPhi_dr.transpose();
  // mGradientPhi_ds_t = mGradientPhi_ds.transpose();
  // mGradientPhi_dt_t = mGradientPhi_dt.transpose();
  
  // save for later
  mInvJac = invJ;
  mInvJacT_x_invJac = invJ.transpose() * invJ;
  mInvJacT = invJ.transpose();
  mDetJac = detJ;
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
  MatrixXd elementStiffnessMatrix(mNumIntPnt,mNumIntPnt);

  mGradientPhi_dx.resize(mNumIntPnt,mNumIntPnt);
  mGradientPhi_dy.resize(mNumIntPnt,mNumIntPnt);
  mGradientPhi_dz.resize(mNumIntPnt,mNumIntPnt);

  mWiDPhi_x.resize(mNumIntPnt,mNumIntPnt);
  mWiDPhi_y.resize(mNumIntPnt,mNumIntPnt);
  mWiDPhi_z.resize(mNumIntPnt,mNumIntPnt);

  // Stiffness Matrix
  // Dx = rx*Dr + sx*Ds + tx*Dt;
  // Dy = ry*Dr + sy*Ds + ty*Dt;
  // Dz = rz*Dr + sz*Ds + tz*Dt;
  // Kxx = detJ*Dx'*Me*Dx
  // Kyy = detJ*Dy'*Me*Dy
  // Kzz = detJ*Dz'*Me*Dz
  // K = Kxx+Kyy+Kzz

  // alternate method
  mGradientPhi_dx =
    (drdx*mGradientPhi_dr +
    dsdx*mGradientPhi_ds +
     dtdx*mGradientPhi_dt).transpose();
  mGradientPhi_dy =
    (drdy*mGradientPhi_dr +
    dsdy*mGradientPhi_ds +
     dtdy*mGradientPhi_dt).transpose();
  mGradientPhi_dz =
    (drdz*mGradientPhi_dr +
    dsdz*mGradientPhi_ds +
     dtdz*mGradientPhi_dt).transpose();
  
  // Me*Dx, Me*Dy, Me*Dz
  RealVec wi_vp2 = mIntegrationWeights * vp2;
  mWiDPhi_x = (wi_vp2).asDiagonal() * mGradientPhi_dx;
  mWiDPhi_y = (wi_vp2).asDiagonal() * mGradientPhi_dy;
  mWiDPhi_z = (wi_vp2).asDiagonal() * mGradientPhi_dz;
  
  // for(int i=0;i<mNumIntPnt;i++) {
  //   mWiDPhi_x.row(i) = mGradientPhi_dx.row(i)*mIntegrationWeights[i];
  //   mWiDPhi_y.row(i) = mGradientPhi_dy.row(i)*mIntegrationWeights[i];
  //   mWiDPhi_z.row(i) = mGradientPhi_dz.row(i)*mIntegrationWeights[i];
  // }
    
  elementStiffnessMatrix = (detJ*mGradientPhi_dx.transpose()*mWiDPhi_x) + (detJ*mGradientPhi_dy.transpose()*mWiDPhi_y) + (detJ*mGradientPhi_dz.transpose()*mWiDPhi_z);
  
  // build transpose as well
  // mGradientPhi_dx_t = mGradientPhi_dx.transpose();
  // mGradientPhi_dy_t = mGradientPhi_dy.transpose();
  // mGradientPhi_dz_t = mGradientPhi_dz.transpose();
  
  return elementStiffnessMatrix;
}



template <typename ConcreteShape>
void Tetrahedra<ConcreteShape>::setupGradientOperator() {

  if(mPlyOrd == 3) {
    mGradientPhi_dr.resize(mNumIntPnt,mNumIntPnt);
    mGradientPhi_ds.resize(mNumIntPnt,mNumIntPnt);
    mGradientPhi_dt.resize(mNumIntPnt,mNumIntPnt);
    dphi_dr_rstn_p3_tetrahedra(mGradientPhi_dr.data());
    dphi_ds_rstn_p3_tetrahedra(mGradientPhi_ds.data());
    dphi_dt_rstn_p3_tetrahedra(mGradientPhi_dt.data());
  } else {        
    std::cerr << "NOT implemented yet!\n";
    MPI::COMM_WORLD.Abort(-1);
  }
    
}

template <typename ConcreteShape>
bool Tetrahedra<ConcreteShape>::attachReceiver(std::unique_ptr<Receiver> &receiver, const bool finalize) {
  printf("TODO: attachedReciever\n");
  exit(1);
}

template <typename ConcreteShape>
double Tetrahedra<ConcreteShape>::CFL_constant() {
  if(mPlyOrd == 3) {
    return 0.375; // by hand (within 10%)    
  } else {
    std::cerr << "ERROR: Order CFL_constant not implemented yet\n";
    exit(1);
  }
}

template <typename ConcreteShape>
double Tetrahedra<ConcreteShape>::estimatedElementRadius() {
  Matrix3d invJ;
  double detJ;

  std::tie(invJ, detJ) = ConcreteShape::inverseJacobian(mVtxCrd);
  Matrix3d J = invJ.inverse();
  VectorXcd eivals = J.eigenvalues();

  // get minimum h (smallest direction)
  Vector3d eivals_norm;
  for(int i=0;i<3;i++) {
    eivals_norm(i) = std::norm(eivals[i]);
  }
  return eivals_norm.minCoeff();
  
}


template <typename ConcreteShape>
VectorXd Tetrahedra<ConcreteShape>::applyTestAndIntegrate(const Ref<const VectorXd> &f) {

  double detJac;
  Matrix3d invJac;
  std::tie(invJac,detJac) = ConcreteShape::inverseJacobian(mVtxCrd);
  return detJac*mIntegrationWeights.array()*f.array();
}

template <typename ConcreteShape>
void Tetrahedra<ConcreteShape>::setBoundaryConditions(std::unique_ptr<Mesh> const &mesh) {
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

/*
 *              (v3)
 *                +--
 *                | \----
 *                |   \--\----
 *                |      \--  \----
 *                |         \-     \---  (v1)
 *                |           \--      \--
 *  ^             |              \--/---  \
 *  |             |             /---\-     \
 *  | (t)         |         /---      \--   \
 *  |             |     /---             \-- \
 *  |      (s)    | /---                    \-\
 *  |       /-    +---------------------------\--
 *  |    /--     (v0)                        (v2)
 *  | /--
 *  +-----------------> (r)
 * 
 * Faces are organized like: bottom (r-s), left (s-t), front (r-t), right (r-s-t)
  * More precisely, faces are composed of vertices and ordered
  * 0: 0 1 2
  * 1: 0 3 1
  * 2: 0 2 3
  * 3: 2 1 3
  * Edge ordering
  * (v0,v1)
  * (v1,v2)
  * (v2,v0)
  * (v0,v3)
  * (v3,v1)
  * (v2,v3)
  */

template <typename ConcreteShape>
std::vector<PetscInt> Tetrahedra<ConcreteShape>::getDofsOnFace(const PetscInt face) {

  std::vector<int> face_ids;
  if(mPlyOrd == 3) {
    switch(face) {
      // see BuildNodesTetrahedraP3() to get surface ids
    case 0: face_ids = {10, 11, 12, 13, 14, 15, 34, 35, 36, 37, 38, 39, 46, 47, 48}; break;
    case 1: face_ids = {16, 17, 18, 19, 20, 21, 34, 35, 40, 41, 42, 43, 46, 47, 49}; break;
    case 3: face_ids = {22, 23, 24, 25, 26, 27, 38, 39, 40, 41, 44, 45, 46, 48, 49}; break;
    case 4: face_ids = {28, 29, 30, 31, 32, 33, 36, 37, 42, 43, 44, 45, 47, 48, 49}; break;
    default: ERROR() << "Only four face on a tetrahedra"; break;
    }
  } else {
    ERROR() << "Not Implemented: getDofsOnFace for Polynomials != 3";
  }
  return face_ids;
}

template <typename ConcreteShape>
std::vector<PetscInt> Tetrahedra<ConcreteShape>::getDofsOnEdge(const PetscInt edge) {

  std::vector<int> edge_ids;
  if(mPlyOrd == 3) {
    switch(edge) {
      // See BuildNodesTetrahedraP3(True) to visualize edge_ids
    case 0: edge_ids = {46,34,35,47}; break;
    case 1: edge_ids = {47,36,37,48}; break;
    case 2: edge_ids = {48,38,39,46}; break;
    case 3: edge_ids = {46,40,41,49}; break;
    case 4: edge_ids = {49,42,43,47}; break;
    case 5: edge_ids = {48,44,45,49}; break;
    default: ERROR() << "Only six edges in a tetrahedra"; break;
    }
  } else {
    ERROR() << "Not Implemented: getDofsOnEdge for Polynomials != 3";
  }
  return edge_ids;
}

template <typename ConcreteShape>
PetscInt Tetrahedra<ConcreteShape>::getDofsOnVtx(const PetscInt vtx) {

  PetscInt vtx_id = -1;
  if(mPlyOrd == 3) {
    switch(vtx) {
    case 0: vtx_id = 46; break;
    case 1: vtx_id = 47; break;
    case 2: vtx_id = 48; break;
    case 3: vtx_id = 49; break;
    default: ERROR() << "Only four vertices in a Tetrahedra"; break;
    }
  } else {
    ERROR() << "Not Implemented: getDofsOnVertex for Polynomials != 3";
  }
  return vtx_id;
}

// Instantiate combinatorical cases.
template class Tetrahedra<TetP1>;
