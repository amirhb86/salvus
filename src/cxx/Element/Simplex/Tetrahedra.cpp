#include <mpi.h>
#include <iostream>
#include <petscdm.h>
#include <petscdmplex.h>
#include <tuple>
#include "Tetrahedra.h"

using namespace Eigen;

// enables std::cout << vector;
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

int Tetrahedra::mNumberVertex = 4;
Eigen::MatrixXd Tetrahedra::mGradientOperator;
Eigen::VectorXd Tetrahedra::mIntegrationWeights;
Eigen::VectorXd Tetrahedra::mIntegrationCoordinates_r;
Eigen::VectorXd Tetrahedra::mIntegrationCoordinates_s;
Eigen::VectorXd Tetrahedra::mIntegrationCoordinates_t;
Eigen::MatrixXd Tetrahedra::mGradientPhi_dr;
Eigen::MatrixXd Tetrahedra::mGradientPhi_ds;
Eigen::MatrixXd Tetrahedra::mGradientPhi_dt;

std::vector<int> getVertsFromPoint(int point, int numVerts, DM &distributed_mesh);

std::tuple<Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd> Tetrahedra::QuadraturePointsForOrder(const int order) {

  if(order == 3) {
    int num_pts = 50;
    Eigen::VectorXd rn(num_pts);
    Eigen::VectorXd sn(num_pts);
    Eigen::VectorXd tn(num_pts);
    coordinates_p3_tetrahedra_rn(rn.data());   
    coordinates_p3_tetrahedra_sn(sn.data());
    coordinates_p3_tetrahedra_tn(tn.data());
    return std::make_tuple(rn,sn,tn);
  }
  else {
    std::cerr << "ERROR: Order NOT implemented!...\n";
    MPI::COMM_WORLD.Abort(-1);
  }
    
}

Eigen::VectorXd Tetrahedra::QuadratureIntegrationWeightForOrder(const int order) {
    
  if(order == 3) {
    int num_pts = 50;
    Eigen::VectorXd wn(num_pts);
    quadrature_weights_p3_tetrahedra(wn.data());
    return wn;
  } else {
    std::cerr << "ERROR: Order NOT implemented!\n";
    MPI::COMM_WORLD.Abort(-1);
  }
}

void tetEdgeHandler(std::vector<std::vector<int>> canonical_edges,
                    std::vector<std::vector<int>> edge_verts,
                    std::vector<int> o,
                 std::vector<std::vector<int>> std_layout,
                    int &dof_counter,
                    VectorXi &element_dof_map
                    );

void tetSurfRemap(int o,
                  int &dof_counter, // current dof count
                  VectorXi &element_dof_map // final dof layout (output)) {
                  );

void tetSurfHandler(std::vector<int> canonical_verts, // the vertices as we expect them
                        std::vector<int> surf_verts, // vertices given from actual surface
                        int o, // mesh given orientation of this surface
                        std::vector<int> std_layout, // expected degrees of freedom layout for this surface
                        int &dof_counter, // current dof count
                        VectorXi &element_dof_map // final dof layout (output)
                    );

VectorXi internalTetMapping(int element, DM &distributed_mesh);

Eigen::VectorXi Tetrahedra::ClosureMapping(const int order, const int dimension,
                                           DM &distributed_mesh) {

  if(order == 3) {

    VectorXi linear_mapping(50);
    for(int i=0;i<50;i++) { linear_mapping[i] = i; }

    auto petsc_mapping = internalTetMapping(mElmNum,distributed_mesh);
    // printf("mapping:\n");
    // for(int i=0;i<25;i++) {
    //   printf("%d:%d\n",i,petsc_mapping[i]);
    // }
    
    VectorXi new_mapping(50);
    for(int i=0;i<50;i++) {
      new_mapping[i] = petsc_mapping[i];
    }
    // return linear_mapping;
    return new_mapping;
    
  } else {
    std::cerr << "ERROR: Order NOT implemented!\n";
    MPI::COMM_WORLD.Abort(-1);
  }
    
}

Eigen::Vector4d Tetrahedra::interpolateAtPoint(double r, double s, double t) {
  // barycentric coordinates l1,l2,l3,l4 (l1+l2+l3+l4=1) give us a linear weighting between coordinates (r,s,t). By computing the barycentric weights, we can compute the interpolated velocity = [l1,l2,l3].dot([v1,v2,v3]) where vN is the vertex velocity.
  // T*[l1;l2;l3] = xyz - v4, where xyz is the desired point, and v4 is the vertex 4
  // thus l123 = T**(-1)*xyz - T**(-1)*v4
  // For reference tet
  // Tref =
  // Matrix([
  // [ 0,  0,  1],
  // [ 0,  1,  0],
  // [-1, -1, -1]])
  // See tetrahedra_element_metrics.py: `interp_vector`
  Eigen::Vector4d interpolator;
  interpolator <<
    -r - s - t + 1, s, r, t;    
  return interpolator;
}

Tetrahedra::Tetrahedra(Options options) {

  mNumDim = 3;

  // Basic properties.
  mPlyOrd = options.PolynomialOrder();
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
        
  // mVtxCrd has 4 vertices
  mElmCtr.resize(3,1);
  mVtxCrd.resize(3,mNumberVertex);
        
  // Integration points and weights
  std::tie(mIntegrationCoordinates_r,mIntegrationCoordinates_s,mIntegrationCoordinates_t) =
    Tetrahedra::QuadraturePointsForOrder(options.PolynomialOrder());        
  mIntegrationWeights = Tetrahedra::QuadratureIntegrationWeightForOrder(options.PolynomialOrder());
          
  setupGradientOperator();
  
}

void Tetrahedra::BuildClosureMapping(DM &distributed_mesh) {

  mClsMap = ClosureMapping(3, 3, distributed_mesh);
  
}

void Tetrahedra::attachVertexCoordinates(DM &distributed_mesh) {

  BuildClosureMapping(distributed_mesh);

  Vec coordinates_local;
  PetscInt coordinate_buffer_size;
  PetscSection coordinate_section;
  PetscReal *coordinates_buffer = NULL;

  DMGetCoordinatesLocal(distributed_mesh, &coordinates_local);
  DMGetCoordinateSection(distributed_mesh, &coordinate_section);
  DMPlexVecGetClosure(distributed_mesh, coordinate_section, coordinates_local, mElmNum,
                      &coordinate_buffer_size, &coordinates_buffer);
  std::vector<PetscReal> coordinates_element(coordinates_buffer, coordinates_buffer+coordinate_buffer_size);
  DMPlexVecRestoreClosure(distributed_mesh, coordinate_section, coordinates_local, mElmNum,
                          &coordinate_buffer_size, &coordinates_buffer);
    
  for (int i = 0; i < mNumberVertex; i++) {        
    mVtxCrd(0,i) = coordinates_element[mNumDim*i+0];
    mVtxCrd(1,i) = coordinates_element[mNumDim*i+1];
    mVtxCrd(2,i) = coordinates_element[mNumDim*i+2];
  }

  // Save element center
  mElmCtr <<
    mVtxCrd.row(0).mean(),
    mVtxCrd.row(1).mean(),
    mVtxCrd.row(2).mean();

}

std::tuple<Eigen::Matrix3d,PetscReal> Tetrahedra::inverseJacobianAtPoint(PetscReal r,
                                                                         PetscReal s,
                                                                         PetscReal t) {

  auto v1x = mVtxCrd(0,0);
  auto v2x = mVtxCrd(0,1);
  auto v3x = mVtxCrd(0,2);
  auto v4x = mVtxCrd(0,3);
  
  auto v1y = mVtxCrd(1,0);
  auto v2y = mVtxCrd(1,1);
  auto v3y = mVtxCrd(1,2);
  auto v4y = mVtxCrd(1,3);
  
  auto v1z = mVtxCrd(2,0);
  auto v2z = mVtxCrd(2,1);
  auto v3z = mVtxCrd(2,2);
  auto v4z = mVtxCrd(2,3);

  // wrong ordering from julia (with vertical normal on bottom)
  // auto dxdr = (v2x-v1x);
  // auto dxds = (v3x-v1x);
  // auto dxdt = (v4x-v1x);
  // auto dydr = (v2y-v1y);
  // auto dyds = (v3y-v1y);
  // auto dydt = (v4y-v1y);
  // auto dzdr = (v2z-v1z);
  // auto dzds = (v3z-v1z);
  // auto dzdt = (v4z-v1z);
  
  // follows correct vertex ordering
  auto dxdr = (v3x-v1x);
  auto dxds = (v2x-v1x);
  auto dxdt = (v4x-v1x);
  auto dydr = (v3y-v1y);
  auto dyds = (v2y-v1y);
  auto dydt = (v4y-v1y);
  auto dzdr = (v3z-v1z);
  auto dzds = (v2z-v1z);
  auto dzdt = (v4z-v1z);
  
  Matrix3d J;
  J <<
    dxdr, dydr, dzdr,
    dxds, dyds, dzds,
    dxdt, dydt, dzdt;

  PetscReal detJ = J.determinant();

  // need to test if this analytical inverse is faster than simply J.inverse();
  
  auto  rx = (1 / detJ) * (dyds*dzdt - dydt*dzds);
  auto  sx = (1 / detJ) * (-dydr*dzdt + dydt*dzdr);
  auto  tx = (1 / detJ) * (dydr*dzds - dyds*dzdr);
  
  auto  ry = (1 / detJ) * (-dxds*dzdt + dxdt*dzds);
  auto  sy = (1 / detJ) * ( dxdr*dzdt - dxdt*dzdr);
  auto  ty = (1 / detJ) * (-dxdr*dzds + dxds*dzdr);
  
  auto  rz = (1 / detJ) * (dxds*dydt - dxdt*dyds);
  auto  sz = (1 / detJ) * (-dxdr*dydt + dxdt*dydr);
  auto  tz = (1 / detJ) * (dxdr*dyds - dxds*dydr);

  Matrix3d inverseJacobian;
  inverseJacobian <<
    rx, sx, tx,
    ry, sy, ty,
    rz, sz, tz;

  return std::make_tuple(inverseJacobian, detJ);
  
}

Eigen::Vector4d Tetrahedra::__attachMaterialProperties(ExodusModel *model, std::string parameter_name) {

  Eigen::Vector4d material_at_vertices(mNumberVertex);

  for (auto i = 0; i < mNumberVertex; i++) {
    material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(
                                                                           mElmCtr, parameter_name, i);
  }
  return material_at_vertices;

}

void Tetrahedra::attachSource(std::vector<std::shared_ptr<Source>> sources) {

  for (auto &source: sources) {
    if (mCheckHull(source->PhysicalLocationX(), source->PhysicalLocationY(), source->PhysicalLocationZ())) {
      Eigen::Vector3d reference_location = inverseCoordinateTransform(source->PhysicalLocationX(),
                                                                      source->PhysicalLocationY(),
                                                                      source->PhysicalLocationZ());
      source->setReferenceLocationR(reference_location(0));
      source->setReferenceLocationS(reference_location(1));
      source->setReferenceLocationT(reference_location(2));
      mSrc.push_back(source);
    }
  }

}

bool Tetrahedra::mCheckHull(double x, double y, double z) {

  auto x1 = mVtxCrd(0,0);
  auto x2 = mVtxCrd(0,1);
  auto x3 = mVtxCrd(0,2);
  auto x4 = mVtxCrd(0,3);
  
  auto y1 = mVtxCrd(1,0);
  auto y2 = mVtxCrd(1,1);
  auto y3 = mVtxCrd(1,2);
  auto y4 = mVtxCrd(1,3);
  
  auto z1 = mVtxCrd(2,0);
  auto z2 = mVtxCrd(2,1);
  auto z3 = mVtxCrd(2,2);
  auto z4 = mVtxCrd(2,3);

  // see tetrahedra_element_metrics.py
  // check barycentric coordinates of the point relative to this
  // triangle. If we violate barycentric coorindate assumptions
  // (l1234 >= 0, l1+l2+l3+l4 = 1) we are not inside the hull.
  double l1 = (x*y2*z3 - x*y2*z4 - x*y3*z2 + x*y3*z4 + x*y4*z2 - x*y4*z3 - x2*y*z3 + x2*y*z4 + x2*y3*z - x2*y3*z4 - x2*y4*z + x2*y4*z3 + x3*y*z2 - x3*y*z4 - x3*y2*z + x3*y2*z4 + x3*y4*z - x3*y4*z2 - x4*y*z2 + x4*y*z3 + x4*y2*z - x4*y2*z3 - x4*y3*z + x4*y3*z2)/(x1*y2*z3 - x1*y2*z4 - x1*y3*z2 + x1*y3*z4 + x1*y4*z2 - x1*y4*z3 - x2*y1*z3 + x2*y1*z4 + x2*y3*z1 - x2*y3*z4 - x2*y4*z1 + x2*y4*z3 + x3*y1*z2 - x3*y1*z4 - x3*y2*z1 + x3*y2*z4 + x3*y4*z1 - x3*y4*z2 - x4*y1*z2 + x4*y1*z3 + x4*y2*z1 - x4*y2*z3 - x4*y3*z1 + x4*y3*z2);
  double l2 = (-x*y1*z3 + x*y1*z4 + x*y3*z1 - x*y3*z4 - x*y4*z1 + x*y4*z3 + x1*y*z3 - x1*y*z4 - x1*y3*z + x1*y3*z4 + x1*y4*z - x1*y4*z3 - x3*y*z1 + x3*y*z4 + x3*y1*z - x3*y1*z4 - x3*y4*z + x3*y4*z1 + x4*y*z1 - x4*y*z3 - x4*y1*z + x4*y1*z3 + x4*y3*z - x4*y3*z1)/(x1*y2*z3 - x1*y2*z4 - x1*y3*z2 + x1*y3*z4 + x1*y4*z2 - x1*y4*z3 - x2*y1*z3 + x2*y1*z4 + x2*y3*z1 - x2*y3*z4 - x2*y4*z1 + x2*y4*z3 + x3*y1*z2 - x3*y1*z4 - x3*y2*z1 + x3*y2*z4 + x3*y4*z1 - x3*y4*z2 - x4*y1*z2 + x4*y1*z3 + x4*y2*z1 - x4*y2*z3 - x4*y3*z1 + x4*y3*z2);
  double l3 = ((x - x4)*((y1 - y4)*((x1 - x4)*(-z2 + z4) + (x2 - x4)*(z1 - z4)) + (z1 - z4)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4))) + (x1 - x4)*(-y*((x1 - x4)*(-z2 + z4) + (x2 - x4)*(z1 - z4)) + y4*((x1 - x4)*(-z2 + z4) + (x2 - x4)*(z1 - z4)) - z*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4)) + z4*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4))))/(-(x1 - x4)*(z3 - z4)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4)) + (x3 - x4)*(z1 - z4)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4)) + ((x1 - x4)*(y3 - y4) - (x3 - x4)*(y1 - y4))*((x1 - x4)*(z2 - z4) - (x2 - x4)*(z1 - z4)));
  double l4 = 1-l1-l2-l3;
    
  if(l1 < 0) return false;
  if(l2 < 0) return false;
  if(l3 < 0) return false;
  if(l4 < 0) return false;
    
  if(fabs(1-(l1+l2+l3)) < 1e-8) return false;

  // no assumptions violated, inside triangle
  return true;
    
}

Eigen::Vector3d Tetrahedra::inverseCoordinateTransform(const double &x, const double &y, const double &z) {

  auto x1 = mVtxCrd(0,0);
  auto x2 = mVtxCrd(0,1);
  auto x3 = mVtxCrd(0,2);
  auto x4 = mVtxCrd(0,3);
  
  auto y1 = mVtxCrd(1,0);
  auto y2 = mVtxCrd(1,1);
  auto y3 = mVtxCrd(1,2);
  auto y4 = mVtxCrd(1,3);
  
  auto z1 = mVtxCrd(2,0);
  auto z2 = mVtxCrd(2,1);
  auto z3 = mVtxCrd(2,2);
  auto z4 = mVtxCrd(2,3);

  // see tetrahedra_element_metrics.py (and https://en.wikipedia.org/wiki/Barycentric_coordinate_system)
  auto r=((x - x4)*(-(y1 - y4)*((x1 - x4)*(z2 - z4) - (x2 - x4)*(z1 - z4)) + (z1 - z4)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4))) - (x1 - x4)*(y*(-(x1 - x4)*(z2 - z4) + (x2 - x4)*(z1 - z4)) - y4*(-(x1 - x4)*(z2 - z4) + (x2 - x4)*(z1 - z4)) + z*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4)) - z4*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4))))/(-(x1 - x4)*(z3 - z4)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4)) + (x3 - x4)*(z1 - z4)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4)) + ((x1 - x4)*(y3 - y4) - (x3 - x4)*(y1 - y4))*((x1 - x4)*(z2 - z4) - (x2 - x4)*(z1 - z4)));

  auto s=(-x*y1*z3 + x*y1*z4 + x*y3*z1 - x*y3*z4 - x*y4*z1 + x*y4*z3 + x1*y*z3 - x1*y*z4 - x1*y3*z + x1*y3*z4 + x1*y4*z - x1*y4*z3 - x3*y*z1 + x3*y*z4 + x3*y1*z - x3*y1*z4 - x3*y4*z + x3*y4*z1 + x4*y*z1 - x4*y*z3 - x4*y1*z + x4*y1*z3 + x4*y3*z - x4*y3*z1)/(x1*y2*z3 - x1*y2*z4 - x1*y3*z2 + x1*y3*z4 + x1*y4*z2 - x1*y4*z3 - x2*y1*z3 + x2*y1*z4 + x2*y3*z1 - x2*y3*z4 - x2*y4*z1 + x2*y4*z3 + x3*y1*z2 - x3*y1*z4 - x3*y2*z1 + x3*y2*z4 + x3*y4*z1 - x3*y4*z2 - x4*y1*z2 + x4*y1*z3 + x4*y2*z1 - x4*y2*z3 - x4*y3*z1 + x4*y3*z2);

  auto t=(-x*y1*z2 + x*y1*z3 + x*y2*z1 - x*y2*z3 - x*y3*z1 + x*y3*z2 + x1*y*z2 - x1*y*z3 - x1*y2*z + x1*y2*z3 + x1*y3*z - x1*y3*z2 - x2*y*z1 + x2*y*z3 + x2*y1*z - x2*y1*z3 - x2*y3*z + x2*y3*z1 + x3*y*z1 - x3*y*z2 - x3*y1*z + x3*y1*z2 + x3*y2*z - x3*y2*z1)/(x1*y2*z3 - x1*y2*z4 - x1*y3*z2 + x1*y3*z4 + x1*y4*z2 - x1*y4*z3 - x2*y1*z3 + x2*y1*z4 + x2*y3*z1 - x2*y3*z4 - x2*y4*z1 + x2*y4*z3 + x3*y1*z2 - x3*y1*z4 - x3*y2*z1 + x3*y2*z4 + x3*y4*z1 - x3*y4*z2 - x4*y1*z2 + x4*y1*z3 + x4*y2*z1 - x4*y2*z3 - x4*y3*z1 + x4*y3*z2);    

  Eigen::Vector3d solution {r, s, t};
  return solution;
}

void Tetrahedra::setupGradientOperator() {

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

// global x-z points on all nodes
std::tuple<Eigen::VectorXd,
           Eigen::VectorXd,
           Eigen::VectorXd> Tetrahedra::buildNodalPoints() {
		
  std::vector<PetscReal> ni(mVtxCrd.size());
	
  Eigen::VectorXd nodalPoints_x(mNumIntPnt);
  Eigen::VectorXd nodalPoints_y(mNumIntPnt);
  Eigen::VectorXd nodalPoints_z(mNumIntPnt);

  auto x1 = mVtxCrd(0,0);
  auto x2 = mVtxCrd(0,1);
  auto x3 = mVtxCrd(0,2);
  auto x4 = mVtxCrd(0,3);
  
  auto y1 = mVtxCrd(1,0);
  auto y2 = mVtxCrd(1,1);
  auto y3 = mVtxCrd(1,2);
  auto y4 = mVtxCrd(1,3);
  
  auto z1 = mVtxCrd(2,0);
  auto z2 = mVtxCrd(2,1);
  auto z3 = mVtxCrd(2,2);
  auto z4 = mVtxCrd(2,3);

  // printf("v0(%f,%f,%f),v1(%f,%f,%f),v2(%f,%f,%f),v3(%f,%f,%f)\n",
  // x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4);
  
  for(auto n = 0; n < mNumIntPnt; n++) {
    auto r = mIntegrationCoordinates_r[n];
    auto s = mIntegrationCoordinates_s[n];
    auto t = mIntegrationCoordinates_t[n];
    // map from reference tet point (r,s,t) to this tet (x,y,z)
    auto xn=r*x3 + s*x2 + t*x4 - x1*(r + s + t - 1);
    auto yn=r*y3 + s*y2 + t*y4 - y1*(r + s + t - 1);
    auto zn=r*z3 + s*z2 + t*z4 - z1*(r + s + t - 1);
    
    nodalPoints_x(n) = xn;
    nodalPoints_y(n) = yn;
    nodalPoints_z(n) = zn;
  }    

  return std::make_tuple(nodalPoints_x,nodalPoints_y,nodalPoints_z);
}

void Tetrahedra::attachReceiver(std::vector<std::shared_ptr<Receiver>> &receivers) {
  printf("TODO: attachedReciever\n");
  exit(1);
}

void Tetrahedra::applyDirichletBoundaries(Mesh *mesh, Options &options, const std::string &fieldname) {

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
