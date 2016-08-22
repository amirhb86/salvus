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
  DMPlexRestoreTransitiveClosure(distributed_mesh, element, PETSC_TRUE, &numPoints, &points);

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

  VectorXi new_mapping;
  if(order == 3) {

    VectorXi linear_mapping(50);
    for(int i=0;i<50;i++) { linear_mapping[i] = i; }

    auto petsc_mapping = internalTetMapping(mElmNum,distributed_mesh);
    // printf("mapping:\n");
    // for(int i=0;i<25;i++) {
    //   printf("%d:%d\n",i,petsc_mapping[i]);
    // }

    new_mapping.setZero(50);
    for(int i=0;i<50;i++) {
      new_mapping[i] = petsc_mapping[i];
    }
    // return linear_mapping;

  } else {
    ERROR() << "Order " << order << " tetrahedra closure mapping not implemented!";
  }
  return new_mapping;
    
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
  
  for(int i=0;i<mNumIntPnt;i++) {

    refGrad.setZero(3);
    for(int j=0;j<mNumIntPnt;j++) {
      refGrad(0) += mGradientPhi_dr(j,i) * field(j);
      refGrad(1) += mGradientPhi_ds(j,i) * field(j);
      refGrad(2) += mGradientPhi_dt(j,i) * field(j);
    }
    
    phyGrad = (mInvJac) * refGrad;
    mGradWork(i,0) = phyGrad(0);
    mGradWork(i,1) = phyGrad(1);
    mGradWork(i,2) = phyGrad(2);
        
  }

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
MatrixXd Tetrahedra<ConcreteShape>::buildStiffnessMatrix(const Ref<const VectorXd>& vp2) {
  
  Eigen::Matrix3d invJ;
  double detJ;
  std::tie(invJ,detJ) = ConcreteShape::inverseJacobian(mVtxCrd);

  mGradientPhi_dr_t = mGradientPhi_dr.transpose();
  mGradientPhi_ds_t = mGradientPhi_ds.transpose();
  mGradientPhi_dt_t = mGradientPhi_dt.transpose();
  
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
        detJ*mIntegrationWeights.dot((vp2.array() * dPhi_dx_i.array() * dPhi_dx_j.array()).matrix()) +
        detJ*mIntegrationWeights.dot((vp2.array() * dPhi_dy_i.array() * dPhi_dy_j.array()).matrix()) +
        detJ*mIntegrationWeights.dot((vp2.array() * dPhi_dz_i.array() * dPhi_dz_j.array()).matrix());

      
      
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
    case 3: edge_ids = {48,38,39,46}; break;
    case 4: edge_ids = {46,40,41,49}; break;
    case 5: edge_ids = {49,42,43,47}; break;
    case 6: edge_ids = {48,44,45,49}; break;
    default: ERROR() << "Only sixe edges in a tetrahedra"; break;
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
    case 3: vtx_id = 48; break;
    case 4: vtx_id = 49; break;
    default: ERROR() << "Only four vertices in a Tetrahedra"; break;
    }
  } else {
    ERROR() << "Not Implemented: getDofsOnVertex for Polynomials != 3";
  }
  return vtx_id;
}


//template <typename ConcreteShape>
//void Tetrahedra<ConcreteShape>::applyDirichletBoundaries(std::unique_ptr<Mesh> const &mesh, std::unique_ptr<Options> const &options,
//                                                         const std::string &fieldname) {
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
template class Tetrahedra<TetP1>;
