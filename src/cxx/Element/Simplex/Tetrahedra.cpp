#include <mpi.h>
#include <iostream>
#include <petscdm.h>
#include <petscdmplex.h>
#include <tuple>
#include "Tetrahedra.h"
#include <Element/Simplex/Tetrahedra/TetP1.h>

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

std::vector<int> getVertsFromPoint(int point, int numVerts, DM &distributed_mesh);

template <typename ConcreteShape>
std::tuple<VectorXd,VectorXd,VectorXd> Tetrahedra<ConcreteShape>::QuadraturePoints(const int order) {

  if(order == 3) {
    int num_pts = 50;
    VectorXd rn(num_pts);
    VectorXd sn(num_pts);
    VectorXd tn(num_pts);
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

template <typename ConcreteShape>
VectorXd Tetrahedra<ConcreteShape>::QuadratureIntegrationWeights(const int order) {
    
  if(order == 3) {
    int num_pts = 50;
    VectorXd wn(num_pts);
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
                 ) {
  
  int num_edges = canonical_edges.size();
    
  for(int i=0;i<num_edges;i++) {
    std::vector<int> oriented_edge_verts(2);
    auto this_edge_verts = edge_verts[i];    
    if(o[i] < 0 ) {
      oriented_edge_verts = {this_edge_verts[0],this_edge_verts[1]};
    }
    else {
      oriented_edge_verts = {this_edge_verts[1],this_edge_verts[0]};      
    }
    std::vector<int> reverse_oriented_edge_verts = {oriented_edge_verts[1],oriented_edge_verts[0]};
    // forward match
    if(oriented_edge_verts == canonical_edges[i]) {
      for(int j=0;j<2;j++) {
        element_dof_map[dof_counter] = std_layout[i][j];
        dof_counter++;
      }
    }
    // reverse match
    else if(canonical_edges[i] == reverse_oriented_edge_verts) {
      element_dof_map[dof_counter] = std_layout[i][1];
      dof_counter++;
      element_dof_map[dof_counter] = std_layout[i][0];
      dof_counter++;      
    }
    else {
      std::cerr << "ERROR: vertices don't match at all, canonical:\n["
                << canonical_edges[i]
                << "], oriented_edge:[\n" << oriented_edge_verts << "]\n";
      exit(1);
    }
  }
}

void tetSurfRemap(int o,
                  int &dof_counter, // current dof count
                  VectorXi &element_dof_map // final dof layout (output)) {
                  ) {
  VectorXi remap(6);
  switch(o) {
  case -3:
    remap << 4, 5, 3, 1, 2, 0;
    break;
  case -2:
    remap << 3, 4, 5, 0, 1, 2;
    break;
  case -1:
    remap << 5, 3, 4, 2, 0, 1;
    break;
  case 0:
    remap << 0, 1, 2, 3, 4, 5;
    break;
  default:
    std::cerr << "Surface orientation not implemented!\n";
    exit(1);    
  }
  for(int i=0;i<6;i++) {
    element_dof_map[dof_counter+i] = dof_counter+remap[i];    
  }
  dof_counter += 6;
  
}

void tetSurfHandler(std::vector<int> canonical_verts, // the vertices as we expect them
                        std::vector<int> surf_verts, // vertices given from actual surface
                        int o, // mesh given orientation of this surface
                        std::vector<int> std_layout, // expected degrees of freedom layout for this surface
                        int &dof_counter, // current dof count
                        VectorXi &element_dof_map // final dof layout (output)
                        ) {

  std::vector<int> oriented_surf_verts;
  if(o < 0) {
    oriented_surf_verts = {surf_verts[2],surf_verts[1],surf_verts[0]};    
  }
  else {
    oriented_surf_verts = {surf_verts[0],surf_verts[1],surf_verts[2]};
  }  
  std::vector<int> map_orig_to_desired(3);
  for(int i=0;i<3;i++) {   
    auto it = std::find(canonical_verts.begin(),canonical_verts.end(),oriented_surf_verts[i]);
    int found_loc = it-canonical_verts.begin();
    // map_orig_to_desired[found_loc] = i;
    map_orig_to_desired[i] = found_loc;
  }
  // std::cout << "surface mapping=" << map_orig_to_desired << "\n";
  // alpha points
  for(int i=0;i<3;i++) {
    element_dof_map[dof_counter] = std_layout[map_orig_to_desired[i]];
    dof_counter++;
  }
  // beta points
  for(int i=0;i<3;i++) {
    element_dof_map[dof_counter] = std_layout[map_orig_to_desired[i]];
    dof_counter++;
  }
}

VectorXi internalTetMapping(int element, DM &distributed_mesh) {

  auto canonical_element_vertices = getVertsFromPoint(element,4,distributed_mesh);
  VectorXi element_dof_map(50);
  element_dof_map.setZero();
  std::vector<std::tuple<int,int>> surfaces(4); //(point_id,orientation)
  std::vector<std::tuple<int,int>> edges(6); //(point_id,orientation)
  int* points = NULL;
  int numPoints;    
  std::vector<std::vector<int>> std_layout = {
    {0,1,2,3,4,5,6,7,8,9},// internal dofs
    // --- surfaces ---
    {10,11,12,13,14,15}, // bottom (r-s)
    {16,17,18,19,20,21}, // "left" (s-t)
    {22,23,24,25,26,27}, // front (r-t)
    {28,29,30,31,32,33}, // "right" (r-s-t)
    // --- edges --- (split into edge components)
    {34,35},{36,37},{38,39}, // bottom (clockwise from inside)
    {40,41},{42,43}, // top t-axis,t-s-slant
    {44,45}, // front r-t-slant
    // --- vertices ---
    {46,47,48,49} // outward normals
  };
  DMPlexGetTransitiveClosure(distributed_mesh,element,PETSC_TRUE,&numPoints,&points);  
  // get surfaces
  for(int i=1;i<5;i++) {
    surfaces[i-1] = std::make_tuple(points[2*i],points[2*i+1]);
  }
  // get edges
  for(int i=5;i<11;i++) {
    edges[i-5] = std::make_tuple(points[2*i],points[2*i+1]);
  }

  int dof_counter = 0;
  // 0-9 are internal points
  for(int i=0;i<10;i++) { element_dof_map[dof_counter] = std_layout[0][i]; dof_counter++; }

  int surf_bottom_id = 1; int surf_left_id = 2; int surf_front_id = 3;
  int surf_right_id = 4;
  int edge_bottom_0 = 5;
  int edge_left_0 = 8;
  int edge_front_0 = 10;
  int vertices_id = 11;

  // assign surface ids
  int id,o;
  // -------------------------
  // ------- bottom ----------
  // -------------------------
  // std:: cout << "computing bottom surface mapping\n";
  
  // std::vector<int> canonical_bottom_verts =
  //   {canonical_element_vertices[0], canonical_element_vertices[1],
  //    canonical_element_vertices[2]};
  std::tie(id,o) = surfaces[0];
  tetSurfRemap(o,dof_counter,element_dof_map);
  // auto bottom_surf_verts = getVertsFromPoint(id,3,distributed_mesh);
  // tetSurfHandler(canonical_bottom_verts,
  //             bottom_surf_verts,
  //             o,
  //             std_layout[surf_bottom_id],
  //             dof_counter,
  //             element_dof_map);
  
  // -------------------------
  // ------- left ----------
  // -------------------------
  // std:: cout << "computing left surface mapping\n";
  // std::vector<int> canonical_left_verts =
  //   {canonical_element_vertices[0], canonical_element_vertices[3],
  //    canonical_element_vertices[1]};
  std::tie(id,o) = surfaces[1];
  tetSurfRemap(o,dof_counter,element_dof_map);
  // auto left_surf_verts = getVertsFromPoint(id,3,distributed_mesh);
  // tetSurfHandler(canonical_left_verts,
  //                left_surf_verts,
  //                o,
  //                std_layout[surf_left_id],
  //                dof_counter,
  //                element_dof_map);
  
  // -------------------------
  // ------- front ----------
  // -------------------------
  // std:: cout << "computing front surface mapping\n";
  // std::vector<int> canonical_front_verts =
  //   {canonical_element_vertices[0], canonical_element_vertices[2],
  //    canonical_element_vertices[3]};
  std::tie(id,o) = surfaces[2];
  tetSurfRemap(o,dof_counter,element_dof_map);
  // auto front_surf_verts = getVertsFromPoint(id,3,distributed_mesh);
  // tetSurfHandler(canonical_front_verts,
  //                front_surf_verts,
  //                o,
  //                std_layout[surf_front_id],
  //                dof_counter,
  //                element_dof_map);
  
  // -------------------------
  // ------- right ----------
  // -------------------------
  // std:: cout << "computing right surface mapping\n";
  // std::vector<int> canonical_right_verts =
  //   {canonical_element_vertices[2], canonical_element_vertices[1],
  //    canonical_element_vertices[3]};
  std::tie(id,o) = surfaces[3];
  tetSurfRemap(o,dof_counter,element_dof_map);
  // auto right_surf_verts = getVertsFromPoint(id,3,distributed_mesh);
  // tetSurfHandler(canonical_right_verts,
  //                right_surf_verts,
  //                o,
  //                std_layout[surf_right_id],
  //                dof_counter,
  //                element_dof_map);
  
  // --------------------------------------------------
  // --------------------- edges ----------------------
  // --------------------------------------------------

  // --------------------------
  // ------ bottom edges ------
  // --------------------------
  // std::cout << "computing bottom edge mapping\n";
  std::vector<std::vector<int>> canonical_bottom_edges =
    {{canonical_element_vertices[0],canonical_element_vertices[1]},
     {canonical_element_vertices[1],canonical_element_vertices[2]},
     {canonical_element_vertices[2],canonical_element_vertices[0]}};               
  std::vector<std::tuple<int,int>> bottom_edge_id_and_orientation =
    {edges[0],edges[1],edges[2]};
  std::vector<int> o_bottom(3);
  std::vector<std::vector<int>> bottom_edge_verts(3);
  std::vector<std::vector<int>> std_layout_edge_bottom(3);
  for(int i=0;i<3;i++) {
    int id,o;
    std::tie(id,o) = bottom_edge_id_and_orientation[i];
    bottom_edge_verts[i] = getVertsFromPoint(id,2,distributed_mesh);
    o_bottom[i] = o;
    std_layout_edge_bottom[i] = std_layout[edge_bottom_0+i];
  }

  tetEdgeHandler(canonical_bottom_edges,
                 bottom_edge_verts,
                 o_bottom,
                 std_layout_edge_bottom,
                 dof_counter,
                 element_dof_map);


  // --------------------------
  // ------ left edges ------
  // --------------------------
  // std::cout << "computing left edge mapping\n";
  std::vector<std::vector<int>> canonical_left_edges =
    {{canonical_element_vertices[0],canonical_element_vertices[3]},
     {canonical_element_vertices[3],canonical_element_vertices[1]}};
  std::vector<std::tuple<int,int>> left_edge_id_and_orientation =
    {edges[3],edges[4]};
  std::vector<int> o_left(2);
  std::vector<std::vector<int>> left_edge_verts(2);
  std::vector<std::vector<int>> std_layout_edge_left(2);
  for(int i=0;i<2;i++) {
    int id,o;
    std::tie(id,o) = left_edge_id_and_orientation[i];
    left_edge_verts[i] = getVertsFromPoint(id,2,distributed_mesh);
    o_left[i] = o;
    std_layout_edge_left[i] = std_layout[edge_left_0+i];
  }
  
  tetEdgeHandler(canonical_left_edges,
                 left_edge_verts,
                 o_left,
                 std_layout_edge_left,
                 dof_counter,
                 element_dof_map);

  // --------------------------
  // ------ front edges ------
  // --------------------------
  // std::cout << "computing front edge mapping\n";
  std::vector<std::vector<int>> canonical_front_edges =
    {{canonical_element_vertices[2],canonical_element_vertices[3]}};
  std::vector<std::tuple<int,int>> front_edge_id_and_orientation =
    {edges[5]};
  std::vector<int> o_front(1);
  std::vector<std::vector<int>> front_edge_verts(1);
  std::vector<std::vector<int>> std_layout_edge_front(1);
  for(int i=0;i<1;i++) {
    int id,o;
    std::tie(id,o) = front_edge_id_and_orientation[i];
    front_edge_verts[i] = getVertsFromPoint(id,2,distributed_mesh);
    o_front[i] = o;
    std_layout_edge_front[i] = std_layout[edge_front_0+i];
  }
  
  tetEdgeHandler(canonical_front_edges,
                 front_edge_verts,
                 o_front,
                 std_layout_edge_front,
                 dof_counter,
                 element_dof_map);

  // vertices
  for(int i=0;i<4;i++) {
    element_dof_map[dof_counter] = std_layout[vertices_id][i];
    dof_counter++;
  }
  return element_dof_map;
}

template <typename ConcreteShape>
VectorXi Tetrahedra<ConcreteShape>::ClosureMapping(const int order, const int dimension,
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

template <typename ConcreteShape>
Tetrahedra<ConcreteShape>::Tetrahedra(Options options) {

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
  
  // Integration points and weights
  std::tie(mIntegrationCoordinates_r,mIntegrationCoordinates_s,mIntegrationCoordinates_t) =
    Tetrahedra<ConcreteShape>::QuadraturePoints(options.PolynomialOrder());
  mIntegrationWeights = Tetrahedra<ConcreteShape>::QuadratureIntegrationWeights(options.PolynomialOrder());
          
  setupGradientOperator();

  mParWork.setZero(mNumIntPnt);
  mStiffWork.setZero(mNumIntPnt);
  mGradWork.setZero(mNumIntPnt, mNumDim);
  
}

template <typename ConcreteShape>
void Tetrahedra<ConcreteShape>::attachVertexCoordinates(DM &distributed_mesh) {

  mClsMap = ClosureMapping(3, 3, distributed_mesh);

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
void Tetrahedra<ConcreteShape>::attachMaterialProperties(const ExodusModel *model, std::string parameter_name) {

  Vector4d material_at_vertices;

  for (auto i = 0; i < mNumVtx; i++) {
    material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(mElmCtr, parameter_name, i);
  }
  
  mPar[parameter_name] = material_at_vertices;

}

template <typename ConcreteShape>
void Tetrahedra<ConcreteShape>::attachSource(std::vector<std::shared_ptr<Source>> sources) {

  for (auto &source: sources) {
    if (ConcreteShape::checkHull(source->PhysicalLocationX(),
                                 source->PhysicalLocationY(),
                                 source->PhysicalLocationZ(),
                                 mVtxCrd)) {
      Vector3d reference_location = ConcreteShape::inverseCoordinateTransform(source->PhysicalLocationX(),
                                                                              source->PhysicalLocationY(),
                                                                              source->PhysicalLocationZ(),
                                                                              mVtxCrd);
      source->setReferenceLocationR(reference_location(0));
      source->setReferenceLocationS(reference_location(1));
      source->setReferenceLocationT(reference_location(2));
      mSrc.push_back(source);
    }
  }

}

template <typename ConcreteShape>
MatrixXd Tetrahedra<ConcreteShape>::computeGradient(const Ref<const VectorXd>& field) {

  // loop version is surprisingly faster.
  // mGradWork.col(0) = mGradientPhi_dx_t*field;
  // mGradWork.col(1) = mGradientPhi_dy_t*field;
  // mGradWork.col(2) = mGradientPhi_dz_t*field;
  
  for(int i=0;i<mNumIntPnt;i++) {
    
    // double grad_x = field.dot(mGradientPhi_dx.col(i));
    // double grad_y = field.dot(mGradientPhi_dy.col(i));
    // double grad_z = field.dot(mGradientPhi_dz.col(i));    
    
    // 15% faster when unrolled
    double grad_x = 0.0;
    double grad_y = 0.0;
    double grad_z = 0.0;

    for(int j=0;j<mNumIntPnt;j++) {
      grad_x += field(j)*mGradientPhi_dx(j,i);
      grad_y += field(j)*mGradientPhi_dy(j,i);
      grad_z += field(j)*mGradientPhi_dz(j,i);
    }
    
    mGradWork(i,0) = grad_x;
    mGradWork(i,1) = grad_y;
    mGradWork(i,2) = grad_z;
    
  }

  return mGradWork;
}

template <typename ConcreteShape>
VectorXd Tetrahedra<ConcreteShape>::applyGradTestAndIntegrate(const Ref<const MatrixXd>& f) {
  
  // auto gfx = f.col(0);
  // auto gfy = f.col(1);
  // auto gfz = f.col(2);
  
  // for(int i=0;i<mNumIntPnt;i++) {

  //   // VectorXd dPhi_dx_i = mGradientPhi_dx.row(i);
  //   // VectorXd dPhi_dy_i = mGradientPhi_dy.row(i);
  //   // VectorXd dPhi_dz_i = mGradientPhi_dz.row(i);
  //   // mStiffWork[i] = detJ*mIntegrationWeights.dot((dPhi_dx_i.array()*gfx.array()
  //   //                                               + dPhi_dy_i.array()*gfy.array()
  //   //                                               + dPhi_dz_i.array()*gfz.array()).matrix());
    
  //   // alternate version using transpose
  //   // mStiffWork[i] = detJ*mIntegrationWeights.dot((mGradientPhi_dx_t.col(i).array()*gfx.array() +
  //   //                                               mGradientPhi_dy_t.col(i).array()*gfy.array() +
  //   //                                               mGradientPhi_dz_t.col(i).array()*gfz.array()).matrix());
    
  //   double ans = 0.0;
  //   for(int j=0;j<mNumIntPnt;j++) {
  //     ans += mIntegrationWeights[i]*(mGradientPhi_dx_t(j,i)*gfx[j] +
  //                                    mGradientPhi_dy_t(j,i)*gfy[j] +
  //                                    mGradientPhi_dz_t(j,i)*gfz[j]);
  //     // ans += mIntegrationWeights[i]*(mGradientPhi_dx(i,j)*f(j,0) +
  //     // mGradientPhi_dy(i,j)*f(j,1) +
  //     // mGradientPhi_dz(i,j)*f(j,2));
  //   }
  //   mStiffWork[i] = mDetJac*ans;
    
  // }  
  
  // return mStiffWork;

  return (mWiDPhi_x*f.col(0) + mWiDPhi_y*f.col(1) + mWiDPhi_z*f.col(2));
  
}

template <typename ConcreteShape>
MatrixXd Tetrahedra<ConcreteShape>::buildStiffnessMatrix(VectorXd velocity) {
  
  Eigen::Matrix3d invJ;
  double detJ;
  std::tie(invJ,detJ) = ConcreteShape::inverseJacobian(mVtxCrd);
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
  
  // loop over matrix(i,j)
  for(int i=0;i<mNumIntPnt;i++) {
      
    Eigen::VectorXd dPhi_dr_i = mGradientPhi_dr.row(i);
    Eigen::VectorXd dPhi_ds_i = mGradientPhi_ds.row(i);
    Eigen::VectorXd dPhi_dt_i = mGradientPhi_dt.row(i);
    auto dPhi_dx_i = dPhi_dr_i*drdx + dPhi_ds_i*dsdx + dPhi_dt_i*dtdx;
    auto dPhi_dy_i = dPhi_dr_i*drdy + dPhi_ds_i*dsdy + dPhi_dt_i*dtdy;
    auto dPhi_dz_i = dPhi_dr_i*drdz + dPhi_ds_i*dsdz + dPhi_dt_i*dtdz;
    mGradientPhi_dx.row(i) = dPhi_dx_i;
    mGradientPhi_dy.row(i) = dPhi_dy_i;
    mGradientPhi_dz.row(i) = dPhi_dz_i;

    mWiDPhi_x.row(i) = detJ * dPhi_dx_i.array() * mIntegrationWeights.array();
    mWiDPhi_y.row(i) = detJ * dPhi_dy_i.array() * mIntegrationWeights.array();
    mWiDPhi_z.row(i) = detJ * dPhi_dz_i.array() * mIntegrationWeights.array();
    
    for(int j=0;j<mNumIntPnt;j++) {
      Eigen::VectorXd dPhi_dr_j = mGradientPhi_dr.row(j);
      Eigen::VectorXd dPhi_ds_j = mGradientPhi_ds.row(j);
      Eigen::VectorXd dPhi_dt_j = mGradientPhi_dt.row(j);
      auto dPhi_dx_j = dPhi_dr_j*drdx + dPhi_ds_j*dsdx + dPhi_dt_j*dtdx;
      auto dPhi_dy_j = dPhi_dr_j*drdy + dPhi_ds_j*dsdy + dPhi_dt_j*dtdy;
      auto dPhi_dz_j = dPhi_dr_j*drdz + dPhi_ds_j*dsdz + dPhi_dt_j*dtdz;
      
      elementStiffnessMatrix(i,j) =
        // with velocity according to model
        detJ*mIntegrationWeights.dot((velocity.array().pow(2) * dPhi_dx_i.array() * dPhi_dx_j.array()).matrix()) +
        detJ*mIntegrationWeights.dot((velocity.array().pow(2) * dPhi_dy_i.array() * dPhi_dy_j.array()).matrix()) +
        detJ*mIntegrationWeights.dot((velocity.array().pow(2) * dPhi_dz_i.array() * dPhi_dz_j.array()).matrix());

      
      
    }
  }

  // build transpose as well
  mGradientPhi_dx_t = mGradientPhi_dx.transpose();
  mGradientPhi_dy_t = mGradientPhi_dy.transpose();
  mGradientPhi_dz_t = mGradientPhi_dz.transpose();
  
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
void Tetrahedra<ConcreteShape>::attachReceiver(std::vector<std::shared_ptr<Receiver>> &receivers) {
  printf("TODO: attachedReciever\n");
  exit(1);
}

template <typename ConcreteShape>
VectorXd Tetrahedra<ConcreteShape>::applyTestAndIntegrate(const Ref<const VectorXd> &f) {

  double detJac;
  Matrix3d invJac;
  std::tie(invJac,detJac) = ConcreteShape::inverseJacobian(mVtxCrd);
  return detJac*mIntegrationWeights.array()*f.array();
}

template <typename ConcreteShape>
void Tetrahedra<ConcreteShape>::setBoundaryConditions(Mesh *mesh) {
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
void Tetrahedra<ConcreteShape>::applyDirichletBoundaries(Mesh *mesh, Options &options, const std::string &fieldname) {

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
template class Tetrahedra<TetP1>;
