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
MatrixXd Hexahedra<ConcreteHex>::mGrd;
template <typename ConcreteHex>
MatrixXd Hexahedra<ConcreteHex>::mGrdT;
template <typename ConcreteHex>
MatrixXd Hexahedra<ConcreteHex>::mGrdWgt;
template <typename ConcreteHex>
MatrixXd Hexahedra<ConcreteHex>::mGrdWgtT;

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

template <typename ConcreteHex>
Hexahedra<ConcreteHex>::Hexahedra(Options options) {

  // Basic properties.
  mPlyOrd = options.PolynomialOrder();
  
  // Gll points.
  mNumDofVtx = 1;
  mNumDofEdg = mPlyOrd - 1;
  mNumDofFac = (mPlyOrd - 1) * (mPlyOrd - 1);
  mNumDofVol = (mPlyOrd - 1) * (mPlyOrd - 1) * (mPlyOrd - 1);

  // Integration points.
  mIntCrdR = Hexahedra<ConcreteHex>::GllPoints(mPlyOrd);
  mIntCrdS = Hexahedra<ConcreteHex>::GllPoints(mPlyOrd);
  mIntCrdT = Hexahedra<ConcreteHex>::GllPoints(mPlyOrd);
  mIntWgtR = Hexahedra<ConcreteHex>::GllIntegrationWeights(mPlyOrd);
  mIntWgtS = Hexahedra<ConcreteHex>::GllIntegrationWeights(mPlyOrd);
  mIntWgtT = Hexahedra<ConcreteHex>::GllIntegrationWeights(mPlyOrd);
  
  // Save number of integration points.
  mNumIntPtsR = mIntCrdR.size();
  mNumIntPtsS = mIntCrdS.size();
  mNumIntPtsT = mIntCrdT.size();
  mNumIntPnt = mNumIntPtsR * mNumIntPtsS * mNumIntPtsT;

  // setup evaluated derivatives of test functions
  mGrd = Hexahedra<ConcreteHex>::setupGradientOperator(mPlyOrd);
  mGrdT = mGrd.transpose();

  mGrdWgt.resize(mNumIntPtsR,mNumIntPtsR);
  for(int i=0;i<mNumIntPtsR;i++) {
    for(int j=0;j<mNumIntPtsR;j++) {
      mGrdWgt(i,j) = mGrd(i,j)*mIntWgtR[i];
    }
  }
  mGrdWgtT = mGrdWgt.transpose();
  
  mDetJac.setZero(mNumIntPnt);
  mParWork.setZero(mNumIntPnt);
  mStiffWork.setZero(mNumIntPnt);
  mGradWork.setZero(mNumIntPnt, mNumDim);

}

// given precomputed point set
std::vector<int> getVertsFromPoint(int point, int numVerts, int numPoints, int* points) {
  std::vector<int> verts(numVerts);
  // vertices are the last `numVerts` entries in points
  for(int i=numPoints-numVerts;i<numPoints;i++) {
    verts[i-(numPoints-numVerts)] = points[2*i];
  }
  return verts;
}

// computes closure every time
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

void edgeHandler(std::vector<std::vector<int>> canonical_edges,
                 std::vector<std::vector<int>> edge_verts,
                 std::vector<int> o,
                 std::vector<std::vector<int>> std_layout,
                 int &dof_counter,
                 VectorXi &element_dof_map
                 ) {
  // Note: bottom & top have 4 edges, front & back only have 2, right
  // & left none because edges are shared between faces.
  int num_edges = canonical_edges.size();
  // std::cout << "canonical_edge_vertices=[";
  // for(int i=0;i<num_edges;i++) {
  //   printf("(%d,%d),",canonical_edges[i][0],canonical_edges[i][1]);
  // }
  // std::cout << "\nthis_edge_vertices=[";
  // for(int i=0;i<num_edges;i++) {
  //   printf("(%d,%d;o=%d),",edge_verts[i][0],edge_verts[i][1],o[i]);
  // }
  // printf("]\n");
  
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
      // for(int j=0;j<2;j++) {
      //   element_dof_map[dof_counter] = std_layout[i][(2-1)-j];
      //   dof_counter++;
      // }
    }
    else {
      std::cerr << "ERROR: vertices don't match at all, canonical:\n["
                << canonical_edges[i]
                << "], oriented_edge:[\n" << oriented_edge_verts << "]\n";
      exit(1);
    }
  }
}

void surfHandler(std::vector<int> canonical_verts, // the vertices as we expect them
                 std::vector<int> surf_verts, // vertices given from actual surface
                 int o, // mesh given orientation of this surface
                 std::vector<int> std_layout, // expected degrees of freedom layout for this surface
                 int &dof_counter, // current dof count
                 VectorXi &element_dof_map // final dof layout (output)
                 ) {

  // std::cout << "canonical_verts=" << canonical_verts << "\n";
  // std::cout << "surf_verts=" << surf_verts << "\n";
  // int this_surf_count = 0;
  std::vector<int> oriented_surf_verts;
  // std::cout << "o=" << o << "\n";
  if(o < 0) {
    oriented_surf_verts = {surf_verts[3],surf_verts[2],surf_verts[1],surf_verts[0]};    
  }
  else {
    oriented_surf_verts = {surf_verts[0],surf_verts[1],surf_verts[2],surf_verts[3]};
  }  
  // std::cout << "oriengted_verts=" << oriented_surf_verts << "\n";
  std::vector<int> map_orig_to_desired(4);
  for(int i=0;i<4;i++) {   
    auto it = std::find(canonical_verts.begin(),canonical_verts.end(),oriented_surf_verts[i]);
    int found_loc = it-canonical_verts.begin();
    // map_orig_to_desired[found_loc] = i;
    map_orig_to_desired[i] = found_loc;
  }
  // std::cout << "mapping:" << map_orig_to_desired << "\n";
  for(int i=0;i<4;i++) {
    element_dof_map[dof_counter] = std_layout[map_orig_to_desired[i]];
    dof_counter++;
  }  
}

// order 3 only
VectorXi internalMappingOrder3(int element, DM &distributed_mesh) {
  
  auto canonical_element_vertices = getVertsFromPoint(element,8,distributed_mesh);
  // std::cout << "canonical_element_vertices=" << canonical_element_vertices << "\n";
  VectorXi element_dof_map(64);
  element_dof_map.setZero();
  std::vector<std::tuple<int,int>> surfaces(6); //(point_id,orientation)
  std::vector<std::tuple<int,int>> edges(12); //(point_id,orientation)
  
  int* points = NULL;
  int numPoints;    
  std::vector<std::vector<int>> std_layout = {
    {0,1,2,3,4,5,6,7},// internal dofs
    // --- surfaces ---
    {8,9,10,11}, // bottom (clockwise from inside)
    {12,13,14,15}, // top (counter-clockwise from outside)
    {16,17,18,19}, // front (counter-clockwise from outside)
    {20,21,22,23}, // back (clockwise from inside)
    {24,25,26,27}, // right (counter-clockwise from outside)
    {28,29,30,31}, // left (clockwise from inside)
    // --- edges --- (split into edge components)
    {32,33},{34,35},{36,37},{38,39}, // bottom (clockwise from inside)
    {40,41},{42,43},{44,45},{46,47}, // top (counter-clockwise from outside)
    {48,49},{50,51}, // front (counter-clockwise from outside)
    {52,53},{54,55}, // back (clockwise from inside)
    // --- vertices ---
    {56,57,58,59,60,61,62,63} // outward normals
  };

  int surf_bottom_id = 1; int surf_top_id = 2; int surf_front_id = 3;
  int surf_back_id = 4; int surf_right_id = 5; int surf_left_id = 6;
  int edge_bottom_0 = 7;
  int edge_top_0 = 11;
  int edge_front_0 = 15;
  int edge_back_0 = 17;
  int vertices_id = 19;
  
  DMPlexGetTransitiveClosure(distributed_mesh,element,PETSC_TRUE,&numPoints,&points);  
  // get surfaces
  for(int i=1;i<7;i++) {
    surfaces[i-1] = std::make_tuple(points[2*i],points[2*i+1]);
  }
  // get edges
  for(int i=7;i<19;i++) {
    edges[i-7] = std::make_tuple(points[2*i],points[2*i+1]);
  }
  
  int dof_counter = 0;
  // 0-7 are internal points
  for(int i=0;i<8;i++) { element_dof_map[dof_counter] = std_layout[0][i]; dof_counter++; }

  // assign surface ids
  int id,o;
  // -------------------------
  // ------- bottom ----------
  // -------------------------
  // std:: cout << "computing bottom surface mapping\n";
  std::vector<int> canonical_bottom_verts =
    {canonical_element_vertices[0], canonical_element_vertices[1],
     canonical_element_vertices[2], canonical_element_vertices[3]};
  std::tie(id,o) = surfaces[0];
  auto bottom_surf_verts = getVertsFromPoint(id,4,distributed_mesh);
  surfHandler(canonical_bottom_verts,
              bottom_surf_verts,
              o,
              std_layout[surf_bottom_id],
              dof_counter,
              element_dof_map);

  // -------------------------
  // ------- top -------------
  // -------------------------
  // std:: cout << "computing top surface mapping\n";
  std::vector<int> canonical_top_verts =
    {canonical_element_vertices[4], canonical_element_vertices[5],
     canonical_element_vertices[6], canonical_element_vertices[7]};
  std::tie(id,o) = surfaces[1];
  auto top_surf_verts = getVertsFromPoint(id,4,distributed_mesh);
  surfHandler(canonical_top_verts,
              top_surf_verts,
              o,
              std_layout[surf_top_id],
              dof_counter,
              element_dof_map);

  // -------------------------
  // ------- front -----------
  // -------------------------
  // std:: cout << "computing front surface mapping\n";
  std::vector<int> canonical_front_verts =
    {canonical_element_vertices[0], canonical_element_vertices[3],
     canonical_element_vertices[5], canonical_element_vertices[4]};
  std::tie(id,o) = surfaces[2];
  auto front_surf_verts = getVertsFromPoint(id,4,distributed_mesh);
  surfHandler(canonical_front_verts,
              front_surf_verts,
              o,
              std_layout[surf_front_id],
              dof_counter,
              element_dof_map);

  // -------------------------
  // ------- back -----------
  // -------------------------
  // std:: cout << "computing back surface mapping\n";
  std::vector<int> canonical_back_verts =
    {canonical_element_vertices[2], canonical_element_vertices[1],
     canonical_element_vertices[7], canonical_element_vertices[6]};
  std::tie(id,o) = surfaces[3];
  auto back_surf_verts = getVertsFromPoint(id,4,distributed_mesh);
  surfHandler(canonical_back_verts,
              back_surf_verts,
              o,
              std_layout[surf_back_id],
              dof_counter,
              element_dof_map);

  // -------------------------
  // ------- right -----------
  // -------------------------
  // std:: cout << "computing right surface mapping\n";
  std::vector<int> canonical_right_verts =
    {canonical_element_vertices[3], canonical_element_vertices[2],
     canonical_element_vertices[6], canonical_element_vertices[5]};
  std::tie(id,o) = surfaces[4];
  auto right_surf_verts = getVertsFromPoint(id,4,distributed_mesh);
  surfHandler(canonical_right_verts,
              right_surf_verts,
              o,
              std_layout[surf_right_id],
              dof_counter,
              element_dof_map);

  // -------------------------
  // ------- left -----------
  // -------------------------
  // std:: cout << "computing left surface mapping\n";
  std::vector<int> canonical_left_verts =
    {canonical_element_vertices[0], canonical_element_vertices[4],
     canonical_element_vertices[7], canonical_element_vertices[1]};
  std::tie(id,o) = surfaces[5];
  auto left_surf_verts = getVertsFromPoint(id,4,distributed_mesh);
  surfHandler(canonical_left_verts,
              left_surf_verts,
              o,
              std_layout[surf_left_id],
              dof_counter,
              element_dof_map);

  // --------------------------
  // ------ bottom edges ------
  // --------------------------
  // std::cout << "computing bottom edge mapping\n";
  std::vector<std::vector<int>> canonical_bottom_edges =
    {{canonical_element_vertices[0],canonical_element_vertices[1]},
     {canonical_element_vertices[1],canonical_element_vertices[2]},
     {canonical_element_vertices[2],canonical_element_vertices[3]},     
     {canonical_element_vertices[3],canonical_element_vertices[0]}};               
  std::vector<std::tuple<int,int>> bottom_edge_id_and_orientation =
    {edges[0],edges[1],edges[2],edges[3]};
  std::vector<int> o_bottom(4);
  std::vector<std::vector<int>> bottom_edge_verts(4);
  std::vector<std::vector<int>> std_layout_edge_bottom(4);
  for(int i=0;i<4;i++) {
    int id,o;
    std::tie(id,o) = bottom_edge_id_and_orientation[i];
    bottom_edge_verts[i] = getVertsFromPoint(id,2,distributed_mesh);
    o_bottom[i] = o;
    std_layout_edge_bottom[i] = std_layout[edge_bottom_0+i];
  }

  edgeHandler(canonical_bottom_edges,
              bottom_edge_verts,
              o_bottom,
              std_layout_edge_bottom,
              dof_counter,
              element_dof_map);

  // --------------------------
  // ------ top edges ---------
  // --------------------------
  // std::cout << "computing top edge mapping\n";
  std::vector<std::vector<int>> canonical_top_edges =
    {{canonical_element_vertices[4],canonical_element_vertices[5]},
     {canonical_element_vertices[5],canonical_element_vertices[6]},
     {canonical_element_vertices[6],canonical_element_vertices[7]},     
     {canonical_element_vertices[7],canonical_element_vertices[4]}};
  std::vector<std::tuple<int,int>> top_edge_id_and_orientation =
    {edges[4],edges[5],edges[6],edges[7]};
  std::vector<int> o_top(4);
  std::vector<std::vector<int>> top_edge_verts(4);
  std::vector<std::vector<int>> std_layout_edge_top(4);
  for(int i=0;i<4;i++) {
    int id,o;
    std::tie(id,o) = top_edge_id_and_orientation[i];
    top_edge_verts[i] = getVertsFromPoint(id,2,distributed_mesh);
    o_top[i] = o;
    std_layout_edge_top[i] = std_layout[edge_top_0+i];
  }
  
  edgeHandler(canonical_top_edges,
              top_edge_verts,
              o_top,
              std_layout_edge_top,
              dof_counter,
              element_dof_map);

  // --------------------------
  // ------ front edges -------
  // --------------------------
  // std::cout << "computing front edge mapping\n";
  std::vector<std::vector<int>> canonical_front_edges =
    {{canonical_element_vertices[3],canonical_element_vertices[5]},
     {canonical_element_vertices[4],canonical_element_vertices[0]}};
  std::vector<std::tuple<int,int>> front_edge_id_and_orientation =
    {edges[8],edges[9]};
  std::vector<int> o_front(2);
  std::vector<std::vector<int>> front_edge_verts(2);
  std::vector<std::vector<int>> std_layout_edge_front(2);
  for(int i=0;i<2;i++) {
    int id,o;
    std::tie(id,o) = front_edge_id_and_orientation[i];
    front_edge_verts[i] = getVertsFromPoint(id,2,distributed_mesh);
    o_front[i] = o;
    std_layout_edge_front[i] = std_layout[edge_front_0+i];
  }
  
  edgeHandler(canonical_front_edges,
              front_edge_verts,
              o_front,
              std_layout_edge_front,
              dof_counter,
              element_dof_map);

  // --------------------------
  // ------ back edges --------
  // --------------------------
  // std::cout << "computing back edge mapping\n";
  std::vector<std::vector<int>> canonical_back_edges =
    {{canonical_element_vertices[1],canonical_element_vertices[7]},
     {canonical_element_vertices[6],canonical_element_vertices[2]}};
  std::vector<std::tuple<int,int>> back_edge_id_and_orientation =
    {edges[10],edges[11]};
  std::vector<int> o_back(2);
  std::vector<std::vector<int>> back_edge_verts(2);
  std::vector<std::vector<int>> std_layout_edge_back(2);
  for(int i=0;i<2;i++) {
    int id,o;
    std::tie(id,o) = back_edge_id_and_orientation[i];
    back_edge_verts[i] = getVertsFromPoint(id,2,distributed_mesh);
    o_back[i] = o;
    std_layout_edge_back[i] = std_layout[edge_back_0+i];
  }
  
  edgeHandler(canonical_back_edges,
              back_edge_verts,
              o_back,
              std_layout_edge_back,
              dof_counter,
              element_dof_map);
  

  // vertices
  for(int i=0;i<8;i++) {
    element_dof_map[dof_counter] = std_layout[vertices_id][i];
    dof_counter++;
  }
  
  // printf("dof_counter=%d\n",dof_counter);
  // std::cout << "element_dof_map=" << element_dof_map << "\n";
  return element_dof_map;
}

template <typename ConcreteHex>
VectorXi Hexahedra<ConcreteHex>::ClosureMapping(int order, int elem_num, DM &distributed_mesh) {
  if( order == 3) {
    auto petsc_mapping = internalMappingOrder3(elem_num,distributed_mesh);
    
    /**
       |  DOFs are ordered: {cell, faces, edges, vertices}
       |  The faces are ordered:
       |   * Bottom
       |   * Top   
       |   * Front 
       |   * Back  
       |   * Right 
       |   * Left  
    */
    
    VectorXi mapping_linear(64);
    for(int i=0;i<64;i++) { mapping_linear[i] = i; }
    
    std::vector<int> interior_mapping =
      {21, 22, 25, 26,
       37, 38, 41, 42};
    std::vector<int> bottom_surface_mapping = {5,9,10,6};
        
    std::vector<int> top_surface_mapping = {53,54,58,57};
    
    std::vector<int> front_surface_mapping = {17,18,34,33};
    
    std::vector<int> back_surface_mapping = {30,29,45,46};

    std::vector<int> right_surface_mapping = {23,27,43,39};

    std::vector<int> left_surface_mapping = {20,36,40,24};
    
    // edges    
    std::vector<int> bottom_edge_mapping = {4,8,13,14,11,7,2,1};
    
    std::vector<int> top_edge_mapping = {49,50,55,59,62,61,56,52};
    
    std::vector<int> front_edge_mapping = {19,35,32,16};
    
    std::vector<int> back_edge_mapping = {28,44,47,31};

    // right and left surface-edges are already covered
    std::vector<int> vertex_mapping = {0,12,15,3,48,51,63,60};

    std::vector<int> full_mapping;
    // concatenate them all together!
    full_mapping.insert(full_mapping.end(),
                        interior_mapping.begin(),
                        interior_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        bottom_surface_mapping.begin(),
                        bottom_surface_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        top_surface_mapping.begin(),
                        top_surface_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        front_surface_mapping.begin(),
                        front_surface_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        back_surface_mapping.begin(),
                        back_surface_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        right_surface_mapping.begin(),
                        right_surface_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        left_surface_mapping.begin(),
                        left_surface_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        bottom_edge_mapping.begin(),
                        bottom_edge_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        top_edge_mapping.begin(),
                        top_edge_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        front_edge_mapping.begin(),
                        front_edge_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        back_edge_mapping.begin(),
                        back_edge_mapping.end());
    full_mapping.insert(full_mapping.end(),
                        vertex_mapping.begin(),
                        vertex_mapping.end());
    Map<VectorXi> mapping(full_mapping.data(),full_mapping.size());    
    
    // MatrixXi viewMap(64,2);
    // for(int i=0;i<64;i++) { viewMap(i,0) = i; }
    // viewMap.col(1) = petsc_mapping;    
    // // std::cout << "Hello!\n";
    // std::cout << "viewMap:\n" << viewMap << "\n";
    // return mapping_linear;
    VectorXi new_mapping(64);
    for(int i=0;i<64;i++) {
      new_mapping[i] = full_mapping[petsc_mapping[i]];
    }
    return new_mapping;
    // return petsc_mapping;
    // return mapping;
    // return mapping_linear;
  }
  
  else {
    std::cerr << "ERROR(ClosureMapping): Order for hexahedra not supported\n";
    VectorXi wrong_mapping(1);
    return wrong_mapping;
    
  }
}

template <typename ConcreteHex>
VectorXd Hexahedra<ConcreteHex>::GllPoints(const int order) {
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

template <typename ConcreteHex>
VectorXd Hexahedra<ConcreteHex>::GllIntegrationWeights(const int order) {
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


template <typename ConcreteHex>
VectorXd Hexahedra<ConcreteHex>::rVectorStride(const Eigen::Ref<const Eigen::VectorXd>& f,
                                                                    const int s_ind, const int t_ind,
                                                                    const int numPtsR, const int numPtsS,
                                                                    const int numPtsT) {

  return Map<const VectorXd>(f.data() + s_ind * numPtsS + t_ind*numPtsS*numPtsR,
                             numPtsR);
}

template <typename ConcreteHex>
VectorXd Hexahedra<ConcreteHex>::sVectorStride(const Eigen::Ref<const Eigen::VectorXd>& f,
                                                                                      const int r_ind, const int t_ind,
                                                                                      const int numPtsR, const int numPtsS,
                                                                                      const int numPtsT) {
  return Map<const VectorXd, 0, InnerStride<>>(f.data() + r_ind + t_ind*numPtsR*numPtsS, numPtsS,
                                               InnerStride<>(numPtsS));
}

template <typename ConcreteHex>
VectorXd Hexahedra<ConcreteHex>::tVectorStride(const Eigen::Ref<const Eigen::VectorXd>& f,
                                                                                      const int r_ind, const int s_ind,
                                                                                      const int numPtsR, const int numPtsS,
                                                                                      const int numPtsT) {
  return Map<const VectorXd, 0, InnerStride<>>(f.data() + r_ind + numPtsR*s_ind, numPtsR,
                                               InnerStride<>(numPtsR*numPtsS));
}

template <typename ConcreteHex>
void Hexahedra<ConcreteHex>::attachVertexCoordinates(Mesh *mesh) {

  // needs building after mesh is loaded
  mClsMap = ClosureMapping(mPlyOrd, mElmNum, mesh->DistributedMesh());
  
  Vec coordinates_local;
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

  for (int i = 0; i < mNumVtx; i++) {
    mVtxCrd(i,0) = coordinates_element[mNumDim * i + 0];
    mVtxCrd(i,1) = coordinates_element[mNumDim * i + 1];
    mVtxCrd(i,2) = coordinates_element[mNumDim * i + 2];
  }

  // Save element center
  mElmCtr <<
    mVtxCrd.col(0).mean(),
    mVtxCrd.col(1).mean(),
    mVtxCrd.col(2).mean();
  
}

template <typename ConcreteHex>
double Hexahedra<ConcreteHex>::CFL_constant() {
  if(mPlyOrd == 3) {
    return 2.4; // determined by hand (about 10% conservative)
  }
  else {
    std::cerr << "ERROR: Order CFL_constant not implemented yet\n";
    exit(1);
  }
}

template <typename ConcreteHex>
double Hexahedra<ConcreteHex>::estimatedElementRadius() {

  Matrix3d invJ;
  double detJ;
  
  Matrix3d invJac;
  Vector3d refGrad;
  int num_pts = mNumIntPtsR*mNumIntPtsS*mNumIntPtsT;
  VectorXd h_pts(num_pts);
  
  // Loop over all GLL points.
  for (int t_ind = 0; t_ind < mNumIntPtsT; t_ind++) {
    for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        // gll index.
        int index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;

        // (r,s,t) coordinates for this point.
        double r = mIntCrdR(r_ind);
        double s = mIntCrdS(s_ind);
        double t = mIntCrdT(t_ind);

        // Optimized gradient for tensorized GLL basis.
        std::tie(invJ, detJ) = ConcreteHex::inverseJacobianAtPoint(r, s, t, mVtxCrd);
        Matrix3d J = invJ.inverse();
        Vector3d eivals_abs = J.eigenvalues().array().abs();
        // get minimum h (smallest direction)

        h_pts(index) = eivals_abs.minCoeff();
      }
    }
  }
  return h_pts.minCoeff();
  
}

template <typename ConcreteHex>
void Hexahedra<ConcreteHex>::attachMaterialProperties(const ExodusModel *model,std::string parameter_name) {

  Eigen::VectorXd material_at_vertices(mNumVtx);

  for (auto i = 0; i < mNumVtx; i++) {
    material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(mElmCtr, parameter_name, i);
  }
  mPar[parameter_name] = material_at_vertices;
  
}

template <typename ConcreteHex>
void Hexahedra<ConcreteHex>::attachReceiver(std::vector<std::shared_ptr<Receiver>> &receivers) {

  for (auto &rec: receivers) {
    double x1 = rec->PysLocX1();
    double x2 = rec->PysLocX2();
    double x3 = rec->PysLocX3();
    if (ConcreteHex::checkHull(x1, x2, x3, mVtxCrd)) {
      Vector3d ref_loc = ConcreteHex::inverseCoordinateTransform(x1, x2, x3, mVtxCrd);
      rec->SetRefLocR(ref_loc(0));
      rec->SetRefLocS(ref_loc(1));
      rec->SetRefLocT(ref_loc(2));
      mRec.push_back(rec);
    }
  }
}

template <typename ConcreteHex>
void Hexahedra<ConcreteHex>::attachSource(std::vector<std::shared_ptr<Source>> sources) {
  for (auto &source: sources) {
    double x1 = source->PhysicalLocationX();
    double x2 = source->PhysicalLocationY();
    double x3 = source->PhysicalLocationZ();
    if (ConcreteHex::checkHull(x1, x2, x3, mVtxCrd)) {
      Vector3d ref_loc = ConcreteHex::inverseCoordinateTransform(x1, x2, x3, mVtxCrd);
      source->setReferenceLocationR(ref_loc(0));
      source->setReferenceLocationS(ref_loc(1));
      source->setReferenceLocationT(ref_loc(2));
      mSrc.push_back(source);
    }
  }
}

template <typename ConcreteHex>
VectorXd Hexahedra<ConcreteHex>::getDeltaFunctionCoefficients(const double r, const double s, const double t) {

  Matrix3d invJ;
  mParWork = interpolateLagrangePolynomials(r, s, t, mPlyOrd);
  for (int t_ind = 0; t_ind < mNumIntPtsT; t_ind++) {
    for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {
      

        double ri = mIntCrdR(r_ind);
        double si = mIntCrdS(s_ind);
        double ti = mIntCrdT(t_ind);

        double detJac;
        std::tie(invJ, detJac) = ConcreteHex::inverseJacobianAtPoint(ri, si, ti, mVtxCrd);
        mParWork(r_ind + s_ind * mNumIntPtsR + t_ind*mNumIntPtsR*mNumIntPtsS) /=
          (mIntWgtR(r_ind) * mIntWgtS(s_ind) * mIntWgtT(t_ind)  * detJac);

      }
    }
  }
  return mParWork;
}

  
template <typename ConcreteHex>
VectorXd Hexahedra<ConcreteHex>::interpolateLagrangePolynomials(const double r,
                                                                   const double s,
                                                                   const double t,
                                                                   const int order)

{

  int n_points = (order + 1) * (order + 1) * (order + 1);
  VectorXd gll_coeffs(n_points);
  if (order == 3) {
    interpolate_order3_hex(r, s, t, gll_coeffs.data());
  }
  else {
    std::cerr << "ERROR: Order not implemented yet\n";
    exit(1);
  }
  
  return gll_coeffs;
  
}

template <typename ConcreteHex>
MatrixXd Hexahedra<ConcreteHex>::setupGradientOperator(const int order) {

  auto rn = Hexahedra<ConcreteHex>::GllPoints(order);
  int num_pts_r = rn.size();
  int num_pts_s = rn.size();
  double r = rn(0);
  double s = rn(0);
  
  MatrixXd grad(num_pts_s, num_pts_r);
  Eigen::MatrixXd test(num_pts_r, num_pts_s);
  for (int i=0; i<num_pts_r; i++) {
    double r = rn[i];
    if (order == 1) {
      interpolate_eps_derivative_order1_square(s, test.data());
    } else if (order == 2) {
      interpolate_eps_derivative_order2_square(r, s, test.data());
    } else if (order == 3) {
      interpolate_eps_derivative_order3_square(r, s, test.data());
    } else if (order == 4) {
      interpolate_eps_derivative_order4_square(r, s, test.data());
    } else if (order == 5) {
      interpolate_eps_derivative_order5_square(r, s, test.data());
    } else if (order == 6) {
      interpolate_eps_derivative_order6_square(r, s, test.data());
    } else if (order == 7) {
      interpolate_eps_derivative_order7_square(r, s, test.data());
    } else if (order == 8) {
      interpolate_eps_derivative_order8_square(r, s, test.data());
    } else if (order == 9) {
      interpolate_eps_derivative_order9_square(r, s, test.data());
    } else if (order == 10) {
      interpolate_eps_derivative_order10_square(r, s, test.data());
    }
    grad.row(i) = test.col(0);
  }
  
  return grad;
}

template <typename ConcreteHex>
MatrixXd Hexahedra<ConcreteHex>::computeGradient(const Ref<const VectorXd> &field) {

  Matrix3d invJac;
  Vector3d refGrad;

  // Loop over all GLL points.
  for (int t_ind = 0; t_ind < mNumIntPtsT; t_ind++) {
    for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        // gll index.
        int index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;

        // (r,s,t) coordinates for this point.
        double r = mIntCrdR(r_ind);
        double s = mIntCrdS(s_ind);
        double t = mIntCrdT(t_ind);

        // Optimized gradient for tensorized GLL basis.
        // std::tie(invJac, mDetJac(index)) = ConcreteHex::inverseJacobianAtPoint(r, s, t, mVtxCrd);
        
        // mGradWork.row(index) = invJac * (refGrad <<
        //                                  mGrd.row(r_ind).dot(rVectorStride(field,s_ind,t_ind,
        //                                                                    mNumIntPtsR,mNumIntPtsS,mNumIntPtsT)),
        //                                  mGrd.row(s_ind).dot(sVectorStride(field, r_ind, t_ind,
        //                                                                    mNumIntPtsR,mNumIntPtsS,mNumIntPtsT)),
        //                                  mGrd.row(t_ind).dot(tVectorStride(field, r_ind, s_ind,
        //                                                                    mNumIntPtsR,mNumIntPtsS,mNumIntPtsT))).finished();
        
        // loop by hand is 10 us faster (27 -> 17)
        refGrad.setZero(3);
        for(int i=0;i<mNumIntPtsR;i++) {
          refGrad(0) += mGrd(r_ind,i)*field(i + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS);
          refGrad(1) += mGrd(s_ind,i)*field(r_ind + i * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS);
          refGrad(2) += mGrd(t_ind,i)*field(r_ind + s_ind * mNumIntPtsR + i * mNumIntPtsR * mNumIntPtsS);
        }
        
        mGradWork.row(index) = mInvJac[index] * refGrad;        


      }
    }
  }

  return mGradWork;

}

template <typename ConcreteHex>
VectorXd Hexahedra<ConcreteHex>::ParAtIntPts(const std::string &par) {

  // Loop over all GLL points.
  for (int t_ind = 0; t_ind < mNumIntPtsT; t_ind++) {
    for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        // gll index.
        int index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;

        // (r,s,t) coordinates for this point.
        double r = mIntCrdR(r_ind);
        double s = mIntCrdS(s_ind);
        double t = mIntCrdT(t_ind);

        mParWork(index) = ConcreteHex::interpolateAtPoint(r,s,t).dot(mPar[par]);
      }
    }
  }

  return mParWork;
}


template <typename ConcreteHex>
VectorXd Hexahedra<ConcreteHex>::applyTestAndIntegrate(const Ref<const VectorXd> &f) {

  int i = 0;
  double detJac;
  Matrix3d invJac;
  VectorXd result(mNumIntPnt);
  for (int t_ind = 0; t_ind < mNumIntPtsS; t_ind++) {
    for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        // gll index.
        int index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;
        
        // (r,s) coordinate at this point.
        double r = mIntCrdR(r_ind);
        double s = mIntCrdS(s_ind);
        double t = mIntCrdT(t_ind);

        std::tie(invJac,detJac) = ConcreteHex::inverseJacobianAtPoint(r,s,t,mVtxCrd);
        result(index) = f(index) * detJac * mIntWgtR(r_ind) * mIntWgtS(s_ind) * mIntWgtS(t_ind);

      }
    }
  }

  return result;

}

template <typename ConcreteHex>
VectorXd Hexahedra<ConcreteHex>::applyGradTestAndIntegrate(const Ref<const MatrixXd> &f) {

  // computes the rotatation into x-y-z, which would normally happen later with more terms.
  Vector3d fi;
  MatrixXd fxyz(f.rows(),3);
  for (int t_ind = 0; t_ind < mNumIntPtsS; t_ind++) {
    for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        // gll index.
        int index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;
        fi << f(index,0),f(index,1),f(index,2);
        fi = mInvJac[index].transpose()*fi;
        fxyz(index,0) = fi[0];
        fxyz(index,1) = fi[1];
        fxyz(index,2) = fi[2];
        
      }
    }
  }

  for (int t_ind = 0; t_ind < mNumIntPtsS; t_ind++) {
    for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {
        
        // map reference gradient (lr,ls,lt) to this element (lx,ly,lz)
        auto lr = mGrd.col(r_ind);
        auto ls = mGrd.col(s_ind);
        auto lt = mGrd.col(t_ind);

        int index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;

        double dphi_r_dfx = 0;
        double dphi_s_dfy = 0;
        double dphi_t_dfz = 0;

        for(int i=0;i<mNumIntPtsR;i++) {
          
          int r_index = i + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;
          int s_index = r_ind + i * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;
          int t_index = r_ind + s_ind * mNumIntPtsR + i * mNumIntPtsR * mNumIntPtsS;
          
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

template <typename ConcreteHex>
void Hexahedra<ConcreteHex>::precomputeConstants() {

  int num_pts = mNumIntPtsR*mNumIntPtsS*mNumIntPtsT;
  mDetJac.resize(num_pts);
  mInvJac.resize(num_pts);
  Matrix3d invJac;
  // Loop over all GLL points.
  for (int t_ind = 0; t_ind < mNumIntPtsT; t_ind++) {
    for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        // gll index.
        int index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;

        // (r,s,t) coordinates for this point.
        double r = mIntCrdR(r_ind);
        double s = mIntCrdS(s_ind);
        double t = mIntCrdT(t_ind);

        // Optimized gradient for tensorized GLL basis.
        std::tie(invJac, mDetJac(index)) = ConcreteHex::inverseJacobianAtPoint(r, s, t, mVtxCrd);
        mInvJac[index] = invJac;
      }
    }
  }
}


template <typename ConcreteHex>
VectorXd Hexahedra<ConcreteHex>::computeStiffnessFull(const Ref<const VectorXd> &field, VectorXd &mVp2) {
  
  Vector3d refGrad;
  Vector3d phyGrad;
  Matrix3d invJac;
  int num_pts = mNumIntPtsR*mNumIntPtsS*mNumIntPtsT;

  VectorXd u_r(num_pts);
  VectorXd u_s(num_pts);
  VectorXd u_t(num_pts);
  MatrixXd f(num_pts,3);
  
  for(int j = 0; j<16; j++) {
    for(int i=0; i<4; i++) {  
      u_r[i+4*j] =
        mGrd(i,0)*field[0+4*j] +
        mGrd(i,1)*field[1+4*j] +
        mGrd(i,2)*field[2+4*j] +
        mGrd(i,3)*field[3+4*j];
    }
  }
  
  // Ds
  for(int k=0;k<4;k++) {
    for(int j=0;j<4;j++) {
      for(int i=0;i<4;i++) {
        u_s[i+4*j+4*4*k] =
          field[i+0*4+4*4*k] * mGrdT(0,j) +
          field[i+1*4+4*4*k] * mGrdT(1,j) +
          field[i+2*4+4*4*k] * mGrdT(2,j) +
          field[i+3*4+4*4*k] * mGrdT(3,j);
      }
    }
  }
  
  // Dt
  for(int j = 0; j<4; j++) {
    for(int i=0; i<16; i++) {  
      u_t[i+16*j] =
        field[i+0*16] * mGrdT(0,j) +
        field[i+1*16] * mGrdT(1,j) +
        field[i+2*16] * mGrdT(2,j) +
        field[i+3*16] * mGrdT(3,j);
      
    }
  }
  
  // Loop over all GLL points.
  for (int t_ind = 0; t_ind < mNumIntPtsT; t_ind++) {
    for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        // gll index.
        int index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;

        // refGrad.setZero(3);
        // for(int i=0;i<mNumIntPtsR;i++) {
        //   refGrad(0) += mGrd(r_ind,i)*field(i + s_ind *
        //                                     mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS);
        //   refGrad(1) += mGrd(s_ind,i)*field(r_ind + i * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS);
        //   refGrad(2) += mGrd(t_ind,i)*field(r_ind + s_ind * mNumIntPtsR + i * mNumIntPtsR * mNumIntPtsS);
        // }
        
        refGrad << u_r[index], u_s[index], u_t[index];
        
        phyGrad = (mInvJac[index].transpose() * mInvJac[index]) * refGrad;
        
        // current version
        f(index,0) = mVp2(index)*mDetJac(index)*phyGrad(0);
        f(index,1) = mVp2(index)*mDetJac(index)*phyGrad(1);
        f(index,2) = mVp2(index)*mDetJac(index)*phyGrad(2);
        
      }
    }
  }

  for(int j = 0; j<16; j++) {
    for(int i=0; i<4; i++) {  
      u_r[i+4*j] =
        mGrdWgtT(i,0)*f(0+4*j,0) +
        mGrdWgtT(i,1)*f(1+4*j,0) +
        mGrdWgtT(i,2)*f(2+4*j,0) +
        mGrdWgtT(i,3)*f(3+4*j,0);
    }
  }
  

  // Ds
  for(int k=0;k<4;k++) {
    for(int j=0;j<4;j++) {
      for(int i=0;i<4;i++) {
        u_s[i+4*j+4*4*k] =
          f(i+0*4+4*4*k,1) * mGrdWgt(0,j) +
          f(i+1*4+4*4*k,1) * mGrdWgt(1,j) +
          f(i+2*4+4*4*k,1) * mGrdWgt(2,j) +
          f(i+3*4+4*4*k,1) * mGrdWgt(3,j);
      }
    }
  }
  
  // Dt
  for(int j = 0; j<4; j++) {
    for(int i=0; i<16; i++) {  
      u_t[i+16*j] =
        f(i+0*16,2) * mGrdWgt(0,j) +
        f(i+1*16,2) * mGrdWgt(1,j) +
        f(i+2*16,2) * mGrdWgt(2,j) +
        f(i+3*16,2) * mGrdWgt(3,j);
      
    }
  }
    
  for (int t_ind = 0; t_ind < mNumIntPtsS; t_ind++) {
    for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {

        int index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;
        // vdev
        mStiffWork(index) = (mIntWgtR(s_ind)*mIntWgtR(t_ind)*u_r[index]+
                             mIntWgtR(r_ind)*mIntWgtR(t_ind)*u_s[index]+
                             mIntWgtR(s_ind)*mIntWgtR(r_ind)*u_t[index]);
        
      }
    }
  }

  return mStiffWork;

  
}


template <typename ConcreteHex>
double Hexahedra<ConcreteHex>::integrateField(const Eigen::Ref<const Eigen::VectorXd> &field) {

  double val = 0;
  Matrix3d inverse_Jacobian;
  double detJ;

  Matrix3d invJac;
  for (int t_ind = 0; t_ind < mNumIntPtsS; t_ind++) {
    for (int s_ind = 0; s_ind < mNumIntPtsS; s_ind++) {
      for (int r_ind = 0; r_ind < mNumIntPtsR; r_ind++) {
        
        // (r,s) coordinate at this point.
        double r = mIntCrdR(r_ind);
        double s = mIntCrdS(s_ind);
        double t = mIntCrdT(t_ind);
        
        std::tie(inverse_Jacobian, detJ) = ConcreteHex::inverseJacobianAtPoint(r, s, t, mVtxCrd);
        int index = r_ind + s_ind * mNumIntPtsR + t_ind * mNumIntPtsR * mNumIntPtsS;
        val += field(index) * mIntWgtR(r_ind) *
          mIntWgtS(s_ind) * mIntWgtT(t_ind) * detJ;
        
      }
    }
  }

  return val;
}

template <typename ConcreteHex>
void Hexahedra<ConcreteHex>::setBoundaryConditions(Mesh *mesh) {
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

template <typename ConcreteHex>
void Hexahedra<ConcreteHex>::applyDirichletBoundaries(Mesh *mesh, Options &options, const std::string &fieldname) {

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


// Instantiate base case.
template class Hexahedra<HexP1>;

