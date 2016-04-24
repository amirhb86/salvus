#include "catch.h"
#include <petsc.h>
#include <Eigen/Dense>
#include <Element/Element.h>
#include <Element/HyperCube/Hexahedra.h>
#include <Element/HyperCube/Hex/AcousticHex.h>

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

std::vector<int> getVertsFromPoint(int point, int numVerts, int numPoints, int* points);
std::vector<int> getVertsFromPoint(int point, int numVerts, Mesh* mesh);

TEST_CASE("Test hexahedra strides","[element/hexahedra]") {
  
  int order = 4;
  Options options;
  options.__SetPolynomialOrder(order);
  AcousticHex ref_hex(options);
  int num_pts = ref_hex.__GetNumIntPtsR() * ref_hex.__GetNumIntPtsS() * ref_hex.__GetNumIntPtsT();
  VectorXd hexDataLayout(num_pts);
  VectorXd hexRn(num_pts);
  VectorXd hexSn(num_pts);
  VectorXd hexTn(num_pts);
  int n=0;
  auto rn = ref_hex.__GetIntCoordR();
  auto sn = ref_hex.__GetIntCoordS();
  auto tn = ref_hex.__GetIntCoordT();
  for(int k=0; k<ref_hex.__GetNumIntPtsT(); k++) {
    double t = tn[k];
    for(int j=0; j<ref_hex.__GetNumIntPtsS(); j++) {
      double s = sn[j];
      for(int i=0; i<ref_hex.__GetNumIntPtsR(); i++) {
        double r = rn[i];
        int index = i + j*ref_hex.__GetNumIntPtsR() + k*ref_hex.__GetNumIntPtsR()*ref_hex.__GetNumIntPtsS();
        hexDataLayout[n] = n;
        hexRn[index] = r;
        hexSn[index] = s;
        hexTn[index] = t;
        n++;
      }
    }
  }
  
  VectorXd index_r = ref_hex.rVectorStride(hexRn,2,4);
  VectorXd index_s = ref_hex.sVectorStride(hexSn,1,3);
  VectorXd index_t = ref_hex.tVectorStride(hexTn,4,2);
  REQUIRE(index_r == ref_hex.__GetIntCoordR());
  REQUIRE(index_s == ref_hex.__GetIntCoordS());
  REQUIRE(index_t == ref_hex.__GetIntCoordT());
  MatrixXi checkIndexST(3,2);
  checkIndexST <<
    0, 0,
    0, 1,
    2, 3;
  
  MatrixXd index_r_n(3,ref_hex.__GetNumIntPtsR());
  index_r_n <<
    0,1,2,3,4,
    25,26,27,28,29,
    85,86,87,88,89;
  MatrixXd index_s_n(3,ref_hex.__GetNumIntPtsS());
  index_s_n <<
    0,5,10,15,20,
    25,30,35,40,45,
    77,82,87,92,97;
  MatrixXd index_t_n(3,ref_hex.__GetNumIntPtsT());
  index_t_n <<
    0,25,50,75,100,
    5,30,55,80,105,
    17,42,67,92,117
    ;

  int r_n, s_n, t_n;
  for(int n=0;n<3;n++) {
    s_n = checkIndexST(n,0);
    t_n = checkIndexST(n,1);
    VectorXd rd = ref_hex.rVectorStride(hexDataLayout,s_n,t_n);
    REQUIRE(rd == index_r_n.row(n).transpose());
    
    r_n = checkIndexST(n,0);
    t_n = checkIndexST(n,1);
    VectorXd sd = ref_hex.sVectorStride(hexDataLayout,r_n,t_n);
    REQUIRE(sd == index_s_n.row(n).transpose());
    
    r_n = checkIndexST(n,0);
    s_n = checkIndexST(n,1);
    VectorXd td = ref_hex.tVectorStride(hexDataLayout,r_n,s_n);
    REQUIRE(td == index_t_n.row(n).transpose());
    
  }
  
}

double five_places(double in) {
  return round(in * 1e2)/1e2;
}

// // given precomputed point set
// std::vector<int> getVertsFromPoint(int point, int numVerts, int numPoints, int* points) {
//   std::vector<int> verts(numVerts);
//   // vertices are the last `numVerts` entries in points
//   for(int i=numPoints-numVerts;i<numPoints;i++) {
//     verts[i-(numPoints-numVerts)] = points[2*i];
//   }
//   return verts;
// }

// // computes closure every time
// std::vector<int> getVertsFromPoint(int point, int numVerts, Mesh* mesh) {
//   int* points = NULL;
//   int numPoints;  
//   DMPlexGetTransitiveClosure(mesh->DistributedMesh(),point,PETSC_TRUE,&numPoints,&points);
//   std::vector<int> verts(numVerts);
//   // vertices are the last `numVerts` entries in points
//   for(int i=numPoints-numVerts;i<numPoints;i++) {
//         verts[i-(numPoints-numVerts)] = points[2*i];
//   }
//   return verts;
// }




// TEST_CASE("Test closure mapping more","[element/hexahedra]") {

//   int order = 3;
//   Options options;
//   options.__SetPolynomialOrder(order);
//   Hexahedra ref_hex(options);
//   int num_pts = ref_hex.__GetNumIntPtsR() * ref_hex.__GetNumIntPtsS() * ref_hex.__GetNumIntPtsT();
//   options.__SetMeshType("newmark");
//   // options.__SetExodusMeshFile("simple_hexmesh_2x1x1_2elements.e");
//   options.__SetExodusMeshFile("simple_hexmesh_2x2x2_8elements.vp4.e");
//   options.__SetElementShape("hex");
//   options.__SetPhysicsSystem("acoustic");

//   // Get mesh.
//   Mesh *mesh = Mesh::factory(options);
//   int mesh_load_option = 1;
//   MatrixXd vertices;
//   MatrixXi cells;
//   if(mesh_load_option == 0) {
//     mesh->read(options);  
//   }
//   else if(mesh_load_option == 1) {
//     vertices.resize(3,12);
//     cells.resize(8,2);
//     MatrixXd vertices_t(12,3);
//     vertices_t <<
//       0, 0, 0,
//       0, 1, 0,
//       1, 1, 0,
//       1, 0, 0,
//       0, 0, 1,
//       1, 0, 1,
//       1, 1, 1,
//       0, 1, 1,
//       2, 0, 0,
//       2, 1, 0,
//       2, 1, 1,
//       2, 0, 1;
//     vertices = vertices_t.transpose();
//     MatrixXi cells_t(2,8);
//     // right-left
//     // cells_t <<
//     //   0, 1, 2, 3, 4, 5, 6, 7,
//     //   3, 2, 9, 8, 5, 11, 10, 6;

//     // vertically stacked (-4)
//     // cells_t <<
//     //   0, 4, 7, 1, 3, 2, 6, 5,
//     //   3, 5, 6, 2, 8, 9, 10, 11;      
    
//     // // vertical with twist (-3)
//     // cells_t <<
//     // 0, 4, 7, 1, 3, 2, 6, 5,
//     // 5, 6, 2, 3, 11, 10, 9, 8;
    
//     // vertical with second twist (-2)
//     // cells_t <<
//     //   0, 4, 7, 1, 3, 2, 6, 5,
//     //   6, 2, 3, 5, 10, 11, 8, 9;
    
//     // vertical with third twist (-1)
//     cells_t <<
//       0, 4, 7, 1, 3, 2, 6, 5,
//       2, 3, 5, 6, 9, 10, 11, 8;
    
//     cells = cells_t.transpose();
//     std::cout << "cells_t=\n" << cells_t << "\n";
      
//     mesh->read(3, 2, 12, 8, cells, vertices);
//     // void Mesh::read(int dim, int numCells, int numVerts, int numVertsPerElem,
//     // Eigen::MatrixXi cells, Eigen::MatrixXd vertex_coords)
//   }
  
//   mesh->setupGlobalDof(ref_hex.NumDofVtx(),
//                        ref_hex.NumDofEdg(),
//                        ref_hex.NumDofFac(),
//                        ref_hex.NumDofVol(),
//                        ref_hex.NumDim());
//   // Register all global field required for time stepping.
//   std::cout << "Field[0]=" << mesh->GlobalFields()[0] << "\n";
//   for (auto field : mesh->GlobalFields()) {        
//     mesh->registerFieldVectors(field);
//   }

//   // Setup reference element.
//   Element *reference_element = Element::factory(options);

//   std::vector<Element*> elements;
//   // Get a list of all local elements.
//   for (int i = 0; i < mesh->NumberElementsLocal(); i++) {
//     elements.push_back(reference_element->clone());
//   }    

//   int numPoints;
//   int* points = NULL;
//   DMPlexGetTransitiveClosure(mesh->DistributedMesh(),0,PETSC_TRUE,&numPoints,&points);

//   auto verts_0 = getVertsFromPoint(0,8,mesh);
//   for(int i=0;i<numPoints;i++) {
//     printf("p[%i]=(%d,%d)\n",i,points[i*2+0],points[i*2+1]);    
//   }

//   for(int i=0;i<8;i++) {
//     printf("vert=%d\n",verts_0[i]);
//   }

//   auto verts_38 = getVertsFromPoint(38,4,mesh);
//   std::cout << "verts surf_38=" << verts_38 << "\n";

//   auto verts_40 = getVertsFromPoint(40,4,mesh);
//   std::cout << "verts surf_40=" << verts_40 << "\n";
  
//   Hexahedra* hex0 = dynamic_cast<Hexahedra*> (elements[0]);
//   // Get vertex coordinates from the PETSc DMPLEX.
//   hex0->SetNum(0);
//   hex0->attachVertexCoordinates(mesh->DistributedMesh());
//   hex0->BuildClosureMapping(mesh);
//   std::cout << "element0 vertices:\n" << hex0->VtxCrd() << "\n";
   
//   // points = NULL;
//   // DMPlexGetTransitiveClosure(mesh->DistributedMesh(),1,PETSC_TRUE,&numPoints,&points);
//   // Hexahedra* hex1 = dynamic_cast<Hexahedra*> (elements[1]);
//   // // Get vertex coordinates from the PETSc DMPLEX.
//   // hex1->SetNum(1);
//   // hex1->attachVertexCoordinates(mesh->DistributedMesh());

//   // for(int i=0;i<numPoints;i++) {
//   //   printf("p[%i]=(%d,%d)\n",i,points[i*2+0],points[i*2+1]);
//   // }
  
//   // std::cout << "element1 vertices:\n" << hex1->VtxCrd() << "\n";
  
//   points = NULL;
//   int Nelem = 1;
  
//   Hexahedra* hexN = dynamic_cast<Hexahedra*> (elements[Nelem]);
//   // Get vertex coordinates from the PETSc DMPLEX.
//   hexN->SetNum(Nelem);
//   hexN->attachVertexCoordinates(mesh->DistributedMesh());
//   hexN->BuildClosureMapping(mesh);

//   DMPlexGetTransitiveClosure(mesh->DistributedMesh(),Nelem,PETSC_TRUE,&numPoints,&points);
//   for(int i=0;i<numPoints;i++) {
//     printf("p[%i]=(%d,%d)\n",i,points[i*2+0],points[i*2+1]);
//   }

//   std::cout << "elementN vertices:\n" << hexN->VtxCrd() << "\n";
  
//   // points = NULL;
//   // int N = 15;
//   // DMPlexGetTransitiveClosure(mesh->DistributedMesh(),N,PETSC_TRUE,&numPoints,&points);

//   // for(int i=0;i<numPoints;i++) {
//   //   printf("p[%i]=(%d,%d)\n",i,points[i*2+0],points[i*2+1]);
//   // }  

  
//   int p_start, p_end;
//   int ierr;
//   for(int d=0;d<5;d++) {
//     ierr = DMPlexGetDepthStratum(mesh->DistributedMesh(), d, &p_start, &p_end);
//     printf("%d=[%d,%d)\n",d,p_start,p_end);
//   }
  
//   int npts = 64;
//   VectorXd one_at_one_node(npts);
//   one_at_one_node.setZero();
//   int selected_point = 15;
//   one_at_one_node[12] = 12;
//   one_at_one_node[13] = 13;
//   one_at_one_node[14] = 14;
//   one_at_one_node[15] = 15;
//   mesh->addFieldFromElement("u", 0,
//                             hex0->ClsMap(),
//                             one_at_one_node);
  
//   one_at_one_node.setZero();
//   one_at_one_node = mesh->getFieldOnElement("u", Nelem,
//                                             hexN->ClsMap());

//   for(int i=0;i<npts;i++) {
//     printf("%d:%f",i,one_at_one_node[i]);
//     if(one_at_one_node[i] == 1) printf("!!!!!! ----  !!!");
//     printf("\n");
//   }
  
  
//   printf("-------------------------------------------\n");
//   printf("-------------------------------------------\n");
//   printf("-------------------------------------------\n");
//   printf("-------------------------------------------\n");
//   printf("-------------------------------------------\n");
//   printf("-------------------------------------------\n");
// }


TEST_CASE("Test closure mapping","[element/hexahedra]") {

  int order = 3;
  
  PetscOptionsClear();
  const char *arg[] = {
    "salvus_test",
    "--testing","true",
    "--exodus_file_name", "simple_hexmesh_2x2x2_8elements.vp4.e",
    "--exodus_model_file_name", "simple_hexmesh_2x2x2_8elements.vp4.e",
    "--mesh_type", "newmark",
    "--element_shape", "hex",
    "--physics_system", "acoustic",
    "--polynomial_order", "3",NULL};
  
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(&argc, &argv, NULL);

  // Set options for exact tests
  Options options;
  options.setOptions();
  
  AcousticHex ref_hex(options);
  int num_pts = ref_hex.__GetNumIntPtsR() * ref_hex.__GetNumIntPtsS() * ref_hex.__GetNumIntPtsT();

  
  // Get mesh.
  Mesh *mesh = Mesh::factory(options);
  mesh->read(options);
  mesh->setupGlobalDof(ref_hex.NumDofVtx(),
                       ref_hex.NumDofEdg(),
                       ref_hex.NumDofFac(),
                       ref_hex.NumDofVol(),
                       ref_hex.NumDim());
  // Register all global field required for time stepping.  
  for (auto field : mesh->GlobalFields()) {        
    mesh->registerFieldVectors(field);
  }
  
  // Setup reference element.
  auto reference_element = Element::factory(options);

  std::vector<std::shared_ptr<Element>> elements;
  // Get a list of all local elements.
  for (int i = 0; i < mesh->NumberElementsLocal(); i++) {
    elements.push_back(std::shared_ptr<Element>(reference_element->clone()));
  }    
  VectorXd expected_multiplier(64);
  expected_multiplier <<
    1, 1, 1, 2, 1, 1, 1, 2,
    1, 1, 1, 2, 2, 2, 2, 4,
    1, 1, 1, 2, 1, 1, 1, 2,
    1, 1, 1, 2, 2, 2, 2, 4,
    1, 1, 1, 2, 1, 1, 1, 2,
    1, 1, 1, 2, 2, 2, 2, 4,
    2, 2, 2, 4, 2, 2, 2, 4,
    2, 2, 2, 4, 4, 4, 4, 8;

  VectorXd value_to_set(64);
  // first 64 prime numbers.
  value_to_set <<
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,
    47,53,59,61,67,71,73,79,83,89,97,101,
    103,107,109,113,127,131,137,139,149,
    151,157,163,167,173,179,181,191,193,
    197,199,211,223,227,229,233,239,241,
    251,257,263,269,271,277,281,283,293,
    307,311;

  std::map<std::tuple<double,double,double>,double> value_lookup;
  
    
  VectorXd computed_value(expected_multiplier.size());
  auto hex0 = std::dynamic_pointer_cast<AcousticHex> (elements[0]);
  // Get vertex coordinates from the PETSc DMPLEX.
  hex0->SetNum(0);
  hex0->attachVertexCoordinates(mesh->DistributedMesh());    
  VectorXd pts0_x,pts0_y,pts0_z;
  std::tie(pts0_x,pts0_y,pts0_z) = hex0->buildNodalPoints();
  hex0->BuildClosureMapping(mesh->DistributedMesh());
  
  for(int i=0;i<pts0_x.size();i++) {
    double x = five_places(pts0_x[i]);
    double y = five_places(pts0_y[i]);
    double z = five_places(pts0_z[i]);
    value_lookup.insert(std::pair<std::tuple<double,double,double>,
                        double>(std::make_tuple(x,y,z),value_to_set[i]));
    // printf("inserted %f @ (%e,%e,%e))\n",value_to_set[i],x,y,z);
  }
  
  for(int n=0;n<expected_multiplier.size();n++) {

    auto x_test = pts0_x[n];
    auto y_test = pts0_y[n];
    auto z_test = pts0_z[n];
    
    VectorXd pts_x,pts_y,pts_z;
    // Set up elements.
    int element_number = 0;
    for (auto &element_gen : elements) {
      auto hex = std::dynamic_pointer_cast<AcousticHex> (element_gen);
      // Give each element a number starting from zero.
      hex->SetNum(element_number++);
      
      // Get vertex coordinates from the PETSc DMPLEX.
      hex->attachVertexCoordinates(mesh->DistributedMesh());
      hex->BuildClosureMapping(mesh->DistributedMesh());
      
      std::tie(pts_x,pts_y,pts_z) = hex->buildNodalPoints();
      // if(hex->Num() == 0) {
      //   std::cout << "pts_x=" << pts_x.transpose() << "\n";
      //   std::cout << "pts_y=" << pts_y.transpose() << "\n";
      //   std::cout << "pts_z=" << pts_z.transpose() << "\n";
      // }
      int desired_nodal_id = -1;
      VectorXd one_at_one_node(pts_x.size());
      one_at_one_node.setZero();
      for(int i=0;i<pts_x.size();i++) {
        // find node that is close;
        if( sqrt(pow(pts_x(i)-x_test,2.0) + pow(pts_y(i) - y_test,2.0) + pow(pts_z(i) - z_test,2.0) ) < 1e-2) {
          desired_nodal_id = i;          
          double x = five_places(pts_x[i]);
          double y = five_places(pts_y[i]);
          double z = five_places(pts_z[i]);
          double value_to_use = value_lookup[std::make_tuple(x,y,z)];
          // printf("Using value: %f\n",value_to_use);
          one_at_one_node[desired_nodal_id] = value_to_use;
        }
      }
      mesh->addFieldFromElement("u", hex->Num(),
                                hex->ClsMap(),
                                one_at_one_node);
    }
    auto hex_test = std::dynamic_pointer_cast<AcousticHex> (elements[0]);
    
    VectorXd one_at_one_node(pts_x.size());
    one_at_one_node.setZero();
    
    one_at_one_node = mesh->getFieldOnElement("u", hex_test->Num(),
                                              hex_test->ClsMap());
    std::tie(pts_x,pts_y,pts_z) = hex_test->buildNodalPoints();
    for(int i=0;i<pts_x.size();i++) {
      // find node that is close;
      if( sqrt(pow(pts_x(i)-x_test,2.0) + pow(pts_y(i) - y_test,2.0) + pow(pts_z(i) - z_test,2.0) ) < 1e-2) {
        computed_value[i] = one_at_one_node[i];
      }
    }
  }
  // std::cout << "expected_value-computed_value=\n" << expected_value-computed_value << "\n";

  VectorXd expected_value = (value_to_set.array() * expected_multiplier.array()).matrix();
    
  REQUIRE((expected_value-computed_value).array().abs().maxCoeff() < 1e-5);
}

TEST_CASE("Test Jacobian mapping","[element/hexahedra]") {

  int order = 3;
  PetscOptionsClear();
  const char *arg[] = {
    "salvus_test",
    "--testing","true",
    "--element_shape", "hex",
    "--physics_system", "acoustic",
    "--polynomial_order", "3",NULL};
  
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(&argc, &argv, NULL);

  // Set options for exact tests
  Options options;
  options.setOptions();
  
  AcousticHex ref_hex(options);
  int num_pts = ref_hex.__GetNumIntPtsR() * ref_hex.__GetNumIntPtsS() * ref_hex.__GetNumIntPtsT();
  auto reference_element = Element::factory(options);  
  auto reference_hex = std::dynamic_pointer_cast<AcousticHex> (reference_element);
  // A distorted element
  Eigen::Matrix<double,3,8> coord;

  // test element starts at (0,0), is 1.0 wide, deep, and tall
  Eigen::Matrix<double,8,3> coord_t;
  coord_t <<
    0,0,0,
    0,1,0,
    1,1,0,
    1,0,0,
    0,0,1,
    1,0,1,
    1,1,1,
    0,1,1;
    
  coord << coord_t.transpose();
  
  reference_hex->SetVtxCrd(coord);

  Matrix3d inv_jacobian_1,inv_jacobian_ref1;
  PetscReal detJ;
  std::tie(inv_jacobian_1,detJ) = reference_hex->inverseJacobianAtPoint(0,0,0);
  inv_jacobian_ref1 <<
    2,0,0,
    0,2,0,
    0,0,2;

  REQUIRE((inv_jacobian_1-inv_jacobian_ref1).array().abs().maxCoeff() < 1e-4);
  REQUIRE((detJ - pow(0.5,3))< 1e-6);

  // test element upside down
  coord_t <<
    0,0,1,
    1,0,1,
    1,1,1,
    0,1,1,
    0,0,0,
    0,1,0,
    1,1,0,
    1,0,0;
    
  coord << coord_t.transpose();
  
  reference_hex->SetVtxCrd(coord);

  Matrix3d inv_jacobian_2,inv_jacobian_ref2;
  std::tie(inv_jacobian_2,detJ) = reference_hex->inverseJacobianAtPoint(0,0,0);
  inv_jacobian_ref2 <<
    0,2,0,
    2,0,0,
    0,0,-2;
  
  REQUIRE((inv_jacobian_2-inv_jacobian_ref2).array().abs().maxCoeff() < 1e-4);
  REQUIRE((detJ - 0.5*0.5*(0.5))< 1e-6);

  
}

TEST_CASE("Test inverse coordinate transform mapping","[element/hexahedra]") {

  int order = 3;

  PetscOptionsClear();
  const char *arg[] = {
    "salvus_test",
    "--testing","true",
    "--element_shape", "hex",
    "--physics_system", "acoustic",
    "--polynomial_order", "3",NULL};
  
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(&argc, &argv, NULL);

  // Set options for exact tests
  Options options;
  options.setOptions();

  auto reference_element = Element::factory(options);  
  auto reference_hex = std::dynamic_pointer_cast<AcousticHex> (reference_element);
  int num_pts = reference_hex->__GetNumIntPtsR() * reference_hex->__GetNumIntPtsS() * reference_hex->__GetNumIntPtsT();  
  // A distorted element
  MatrixXd coord(3,8);

  // test element starts at (0,0), is 1.0 wide and deep, 0.5 tall, and has top face shifted by (+1/2,+1/2,0) on all vertices, so that it is slightly slanted.
  coord <<
    0.0, 0.0, 1.0, 1.0, 1.0/2.0, 3.0/2.0, 3.0/2.0, 1.0/2.0, // x
    0.0, 1.0, 1.0, 0.0, 1.0/2.0, 1.0/2.0, 3.0/2.0, 3.0/2.0, // y
    0.0, 0.0, 0.0, 0.0 ,1.0/2.0, 1.0/2.0, 1.0/2.0, 1.0/2.0; // z

  reference_hex->SetVtxCrd(coord);
  VectorXd pts_x,pts_y,pts_z;
  std::tie(pts_x,pts_y,pts_z) = reference_hex->buildNodalPoints();

  VectorXd pts_r(num_pts);
  VectorXd pts_s(num_pts);
  VectorXd pts_t(num_pts);

  VectorXd ref_pts_r(num_pts);
  VectorXd ref_pts_s(num_pts);
  VectorXd ref_pts_t(num_pts);
  
  int idx = 0;
  for (auto k = 0; k < reference_hex->__GetNumIntPtsT(); k++) {
    for (auto j = 0; j < reference_hex->__GetNumIntPtsS(); j++) {      
        for (auto i = 0; i < reference_hex->__GetNumIntPtsR(); i++) {
        
        double r = reference_hex->__GetIntCoordR()[i];
        double s = reference_hex->__GetIntCoordS()[j];
        double t = reference_hex->__GetIntCoordT()[k];
        ref_pts_r[idx] = r;
        ref_pts_s[idx] = s;
        ref_pts_t[idx] = t;
        
        auto vec3 = reference_hex->inverseCoordinateTransform(pts_x[idx],pts_y[idx],pts_z[idx]);
        pts_r[idx] = vec3[0];
        pts_s[idx] = vec3[1];
        pts_t[idx] = vec3[2];

        // printf("(%f,%f,%f) =?= (%f,%f,%f)\n",r,s,t,pts_r[idx],pts_s[idx],pts_t[idx]);
        
        idx++;
      }
    }
  }
  
  REQUIRE((pts_r-ref_pts_r).array().abs().maxCoeff() < 1e-4);
  REQUIRE((pts_s-ref_pts_s).array().abs().maxCoeff() < 1e-4);
  REQUIRE((pts_t-ref_pts_t).array().abs().maxCoeff() < 1e-4);
  
}

TEST_CASE("Test reference mapping","[element/hexahedra]") {

  PetscOptionsClear();
  const char *arg[] = {
    "salvus_test",
    "--testing","true",
    "--element_shape", "hex",
    "--physics_system", "acoustic",
    "--polynomial_order", "3",NULL};
  int order = 3;
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(&argc, &argv, NULL);
  
  // Set options for exact tests
  Options options;
  options.setOptions();
  
  AcousticHex ref_hex(options);
  int num_pts = ref_hex.__GetNumIntPtsR() * ref_hex.__GetNumIntPtsS() * ref_hex.__GetNumIntPtsT();  
  auto reference_element = Element::factory(options);  
  auto reference_hex = std::dynamic_pointer_cast<AcousticHex> (reference_element);
  // A distorted element
  Eigen::Matrix<double,3,8> coord;

  // test element starts at (0,0), is 1.0 wide and deep, 0.5 tall, and has top face shifted by (+1/2,+1/2,0) on all vertices, so that it is slightly slanted.
  coord <<
    0,0,1,1,1.0/2.0,3.0/2.0,3.0/2.0,1.0/2.0, // x
    0,1,1,0,1.0/2.0,1.0/2.0,3.0/2.0,3.0/2.0, // y
    0,0,0,0,1.0/2.0,1.0/2.0,1.0/2.0,1.0/2.0; // z
  
  reference_hex->SetVtxCrd(coord);

  VectorXd pts_x,pts_y,pts_z;
  std::tie(pts_x,pts_y,pts_z) = reference_hex->buildNodalPoints();

  VectorXd ref_pts_x(num_pts);
  VectorXd ref_pts_y(num_pts);
  VectorXd ref_pts_z(num_pts);
  
  int idx = 0;
  for (auto k = 0; k < ref_hex.__GetNumIntPtsT(); k++) {
    for (auto j = 0; j < ref_hex.__GetNumIntPtsS(); j++) {      
        for (auto i = 0; i < ref_hex.__GetNumIntPtsR(); i++) {
        
        double r = ref_hex.__GetIntCoordR()[i];
        double s = ref_hex.__GetIntCoordS()[j];
        double t = ref_hex.__GetIntCoordT()[k];

        ref_pts_x[idx] = (r+1)/2.0 + 0.5*(t+1)/2.0;
        ref_pts_y[idx] = (s+1)/2.0 + 0.5*(t+1)/2.0;
        ref_pts_z[idx] = (t+1)/4.0;
        
        idx++;
      }
    }
  }
  
  REQUIRE((pts_x-ref_pts_x).array().abs().maxCoeff() < 1e-4);
  REQUIRE((pts_y-ref_pts_y).array().abs().maxCoeff() < 1e-4);
  REQUIRE((pts_z-ref_pts_z).array().abs().maxCoeff() < 1e-4);
  
}

TEST_CASE("Test hex velocity interpolation", "[element]") {

  Eigen::MatrixXd mMaterialVelocityAtVertices_i(3,8);
  mMaterialVelocityAtVertices_i <<
    1,1,2,2,1,2,2,1, // left right
    1,2,2,1,1,1,2,2, // front back
    1,1,1,1,2,2,2,2; // bottom top

    
  Eigen::MatrixXd check_velocity_i(3,27);
  check_velocity_i <<
    1.0,1.5,2.0,1.0,1.5,2.0,1.0,1.5,2.0,// left right lowest
    1.0,1.5,2.0,1.0,1.5,2.0,1.0,1.5,2.0,// left right middle
    1.0,1.5,2.0,1.0,1.5,2.0,1.0,1.5,2.0,// left right highest
    // ---
    1.0,1.0,1.0,1.5,1.5,1.5,2.0,2.0,2.0,// front back lowest 
    1.0,1.0,1.0,1.5,1.5,1.5,2.0,2.0,2.0,// front back middle 
    1.0,1.0,1.0,1.5,1.5,1.5,2.0,2.0,2.0,// front back highest
    // ---
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,// bottom top lowest 
    1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,// bottom top middle 
    2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0;// bottom top highest
    
    
  for(int i=0;i<mMaterialVelocityAtVertices_i.rows();i++) {
    Eigen::VectorXd mMaterialVelocityAtVertices = mMaterialVelocityAtVertices_i.row(i);
    Eigen::VectorXd check_velocity = check_velocity_i.row(i);
    Options options;
    options.__SetPolynomialOrder(2);
    
    auto mIntegrationCoordinatesR = Hexahedra::GllPointsForOrder(options.PolynomialOrder());
    auto mIntegrationCoordinatesS = Hexahedra::GllPointsForOrder(options.PolynomialOrder());
    auto mIntegrationCoordinatesT = Hexahedra::GllPointsForOrder(options.PolynomialOrder());
    auto numpts_r = mIntegrationCoordinatesR.size();
    auto numpts_s = mIntegrationCoordinatesS.size();
    auto numpts_t = mIntegrationCoordinatesT.size();
    auto numpts = numpts_r*numpts_s*numpts_t;
    
    int n=0;
    Eigen::VectorXd velocity(numpts);
    for (auto t_index = 0; t_index < numpts_t; t_index++) {
      for (auto s_index = 0; s_index < numpts_s; s_index++) {
        for (auto r_index = 0; r_index < numpts_r; r_index++) {

          // R and s coordinates.
          double t = mIntegrationCoordinatesT[t_index];
          double s = mIntegrationCoordinatesS[s_index];
          double r = mIntegrationCoordinatesR[r_index];
          // Get material parameters at this node.
          // new way
          auto interpolate1 = Hexahedra::interpolateAtPoint(r, s, t);              
          velocity[n] = interpolate1.dot(mMaterialVelocityAtVertices);
          n++;
        }
      }
    }
    
    REQUIRE((velocity.array()-check_velocity.array()).abs().maxCoeff() < 1e-5);
  }
}
