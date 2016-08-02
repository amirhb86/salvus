#include <mpi.h>
#include <Mesh/Mesh.h>
#include <Source/Source.h>
#include <Utilities/Options.h>
#include <Model/ExodusModel.h>
#include <Element/Simplex/TriP1.h>
#include <Element/Simplex/Triangle.h>

extern "C" {
#include <Element/Simplex/Autogen/p3_triangle.h>
}

using namespace Eigen;
/*
 * STATIC variables WHICH ARE ONLY ON THE REFERENCE ELEMENT.
 */
template <typename ConcreteShape>
MatrixXd Triangle<ConcreteShape>::mGradientPhi_dr;
template <typename ConcreteShape>
MatrixXd Triangle<ConcreteShape>::mGradientPhi_ds;

template <typename ConcreteShape>    
std::tuple<VectorXd,VectorXd> Triangle<ConcreteShape>::QuadraturePoints(const int order) {

  VectorXd rn, sn;
  if(order == 3) {
    int num_pts = 12;
    rn.setZero(num_pts);
    sn.setZero(num_pts);
    coordinates_p3_triangle_rn(rn.data());
    coordinates_p3_triangle_sn(sn.data());
  }
  else {
    std::cerr << "ERROR: Order NOT implemented!...\n";
    MPI::COMM_WORLD.Abort(-1);
  }
  return std::make_tuple(rn,sn);
}

template <typename ConcreteShape>
VectorXd Triangle<ConcreteShape>::QuadratureIntegrationWeight(const int order) {

  VectorXd wn;
  if(order == 3) {
    int num_pts = 12;
    wn.setZero(num_pts);
    quadrature_weights_p3_triangle(wn.data());
  } else {
    std::cerr << "ERROR: Order NOT implemented!\n";
    MPI::COMM_WORLD.Abort(-1);
  }
  return wn;
}

template <typename ConcreteShape>
VectorXi Triangle<ConcreteShape>::ClosureMapping(const int order, const int dimension) {

  VectorXi closure;
  if(order == 3) {

    int num_pts = 12;

    closure.setZero(num_pts);
    closure << 0,1,2,3,4,5,6,7,8,9,10,11;

  } else {
    std::cerr << "ERROR: Order NOT implemented!\n";
    MPI::COMM_WORLD.Abort(-1);
  }
  return closure;
    
}

template <typename ConcreteShape>
Triangle<ConcreteShape>::Triangle(std::unique_ptr<Options> const &options) {

  // Basic properties.
  mPlyOrd = options->PolynomialOrder();
  if(mPlyOrd == 3) {
    mNumIntPnt = 12;
    mNumDofEdg = 2;
    mNumDofFac = 3;
    mNumDofVol = 0;
  }
  else {
    std::cerr << "ERROR: Order NOT implemented!\n";
    MPI::COMM_WORLD.Abort(-1);
  }
        
  // Nodal collocation points on edge and surface
  mNumDofVtx = 1;
        
  // Integration points and weights
  std::tie(mIntegrationCoordinates_r,mIntegrationCoordinates_s) =
    Triangle<ConcreteShape>::QuadraturePoints(options->PolynomialOrder());
  mIntegrationWeights = Triangle<ConcreteShape>::QuadratureIntegrationWeight(options->PolynomialOrder());
        
  mClsMap = Triangle<ConcreteShape>::ClosureMapping(options->PolynomialOrder(), mNumDim);
  setupGradientOperator();

  mDetJac.setZero(mNumIntPnt);
  mParWork.setZero(mNumIntPnt);
  mStiffWork.setZero(mNumIntPnt);
  mGradWork.setZero(mNumIntPnt, mNumDim);
  
}

template <typename ConcreteShape>
void Triangle<ConcreteShape>::attachVertexCoordinates(std::unique_ptr<Mesh> const &mesh) {

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
  }

  // Save element center
  mElmCtr <<
    mVtxCrd.col(0).mean(),
    mVtxCrd.col(1).mean();
}

template <typename ConcreteShape>
double Triangle<ConcreteShape>::CFL_constant() {
  if(mPlyOrd == 3) {
    return 4.5; // by hand (within 10%)
  } else {
    std::cerr << "ERROR: Order CFL_constant not implemented yet\n";
    exit(1);
  }
}

template <typename ConcreteShape>
double Triangle<ConcreteShape>::estimatedElementRadius() {
  Matrix2d invJ;
  double detJ;

  std::tie(invJ, detJ) = ConcreteShape::inverseJacobian(mVtxCrd);
  Matrix2d J = invJ.inverse();
  VectorXcd eivals = J.eigenvalues();

  // get minimum h (smallest direction)
  Vector2d eivals_norm;
  for(int i=0;i<2;i++) {
    eivals_norm(i) = std::norm(eivals[i]);
  }
  return eivals_norm.minCoeff();
  
}


template <typename ConcreteShape>
void Triangle<ConcreteShape>::attachMaterialProperties(std::unique_ptr<ExodusModel> const &model, std::string parameter_name) {

  Vector3d material_at_vertices;

  for (auto i = 0; i < mNumVtx; i++) {
    material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(mElmCtr, parameter_name, i);
  }
  mPar[parameter_name] = material_at_vertices;
}

template <typename ConcreteShape>
bool Triangle<ConcreteShape>::attachSource(std::unique_ptr<Source> &source, const bool finalize) {

  if(!source) { return false; }
  if (ConcreteShape::checkHull(source->LocX(), source->LocZ(), mVtxCrd)) {
    if (!finalize) { return true; }
    Vector2d reference_location = ConcreteShape::inverseCoordinateTransform(source->LocX(),
                                                                            source->LocZ(),
                                                                            mVtxCrd);
    source->SetLocR(reference_location(0));
    source->SetLocS(reference_location(1));
    mSrc.push_back(std::move(source));
    return true;
  }
  return false;
}

template <typename ConcreteShape>
MatrixXd Triangle<ConcreteShape>::buildStiffnessMatrix(VectorXd velocity) {

  Matrix2d invJ;
  double detJ;
  std::tie(invJ,detJ) = ConcreteShape::inverseJacobian(mVtxCrd);
  //Jinv= rx, sx,
  //      rz, sz;
  auto drdx = invJ(0,0);
  auto dsdx = invJ(0,1);
  auto drdz = invJ(1,0);
  auto dsdz = invJ(1,1);

  // build material on all nodes
  MatrixXd elementStiffnessMatrix(mNumIntPnt,mNumIntPnt);
    
  // interpolate velocity at all nodes
  // for(int i=0;i<mNumIntPnt;i++) {
  //   auto r = mIntegrationCoordinates_r[i];
  //   auto s = mIntegrationCoordinates_s[i];
  //   velocity(i) = interpolateAtPoint(r,s).dot(mMaterialVelocityAtVertices);
  // }  
  
  // loop over matrix(i,j)
  for(int i=0;i<mNumIntPnt;i++) {
      
    VectorXd dPhi_dr_i = mGradientPhi_dr.row(i);
    VectorXd dPhi_ds_i = mGradientPhi_ds.row(i);
    auto dPhi_dx_i = dPhi_dr_i*drdx + dPhi_ds_i*dsdx;
    auto dPhi_dz_i = dPhi_dr_i*drdz + dPhi_ds_i*dsdz;
    for(int j=0;j<mNumIntPnt;j++) {
      VectorXd dPhi_dr_j = mGradientPhi_dr.row(j);
      VectorXd dPhi_ds_j = mGradientPhi_ds.row(j);
      auto dPhi_dx_j = dPhi_dr_j*drdx + dPhi_ds_j*dsdx;
      auto dPhi_dz_j = dPhi_dr_j*drdz + dPhi_ds_j*dsdz;
            
      elementStiffnessMatrix(i,j) = detJ*mIntegrationWeights.dot((velocity.array().pow(2) * dPhi_dx_i.array() * dPhi_dx_j.array()).matrix()) +
        detJ*mIntegrationWeights.dot((velocity.array().pow(2) * dPhi_dz_i.array() * dPhi_dz_j.array()).matrix());
    }
  }
  return elementStiffnessMatrix;
}


template <typename ConcreteShape>
void Triangle<ConcreteShape>::setupGradientOperator() {

  if(mPlyOrd == 3) {
    mGradientPhi_dr.resize(mNumIntPnt,mNumIntPnt);
    mGradientPhi_ds.resize(mNumIntPnt,mNumIntPnt);
    dphi_dr_rsn_p3_triangle(mGradientPhi_dr.data());
    dphi_ds_rsn_p3_triangle(mGradientPhi_ds.data());
  } else {        
    std::cerr << "NOT implemented yet!\n";
    MPI::COMM_WORLD.Abort(-1);
  }
}

template <typename ConcreteShape>
double Triangle<ConcreteShape>::integrateField(const VectorXd &field) {

  double val = 0;
  Matrix<double,2,2> inverse_Jacobian;
  double detJ;
  std::tie(inverse_Jacobian,detJ) = ConcreteShape::inverseJacobian(mVtxCrd);
  val = detJ*field.dot(mIntegrationWeights);
  return val;

}

template <typename ConcreteShape>
bool Triangle<ConcreteShape>::attachReceiver(std::unique_ptr<Receiver> &receiver, const bool finalize) {
  printf("TODO: attachedReciever\n");
  exit(1);
}

template <typename ConcreteShape>
VectorXd Triangle<ConcreteShape>::applyTestAndIntegrate(const Ref<const VectorXd> &f) {

  double detJac;
  Matrix2d invJac;
  std::tie(invJac,detJac) = ConcreteShape::inverseJacobian(mVtxCrd);
  return detJac*mIntegrationWeights.array()*f.array();
}

template <typename ConcreteShape>
void Triangle<ConcreteShape>::setBoundaryConditions(std::unique_ptr<Mesh> const &mesh) {
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
VectorXd Triangle<ConcreteShape>::getDeltaFunctionCoefficients(const double r, const double s) {
  std::cerr << "ERROR: Not implemented getDeltaFunctionCoefficients\n";
  exit(1);
}

template <typename ConcreteShape>
VectorXd Triangle<ConcreteShape>::ParAtIntPts(const std::string &par) {

  // interpolate velocity at all nodes
  for(int i=0;i<mNumIntPnt;i++) {
    auto r = mIntegrationCoordinates_r[i];
    auto s = mIntegrationCoordinates_s[i];
    mParWork(i) = ConcreteShape::interpolateAtPoint(r,s).dot(mPar[par]);
  }  
  
  // for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
  //   for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

  //     double r = mIntCrdR(r_ind);
  //     double s = mIntCrdS(s_ind);
  //     mParWork(r_ind + s_ind*mNumIntPtsR) = ConcreteShape::interpolateAtPoint(r,s).dot(mPar[par]);

  //   }
  // }
  
  return mParWork;
}


//template <typename ConcreteShape>
//void Triangle<ConcreteShape>::applyDirichletBoundaries(std::unique_ptr<Mesh> const &mesh, std::unique_ptr<Options> const &options,
//                                                       const std::string &fieldname) {
//
//  if (! mBndElm) return;
//
//  double value = 0;
//  auto dirchlet_boundary_names = options->DirichletBoundaries();
//  for (auto &bndry: dirchlet_boundary_names) {
//    auto faceids = mBnd[bndry];
//    for (auto &faceid: faceids) {
//      auto field = mesh->getFieldOnFace(fieldname, faceid);
//      field = 0 * field.array() + value;
//      mesh->setFieldFromFace(fieldname, faceid, field);
//    }
//  }
//}

// Instantiate combinatorical cases.
template class Triangle<TriP1>;
