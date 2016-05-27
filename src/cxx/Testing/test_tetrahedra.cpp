#include "catch.h"
#include <petsc.h>
#include <Eigen/Dense>
#include <Element/Simplex/Tetrahedra.h>
#include <Element/Simplex/TetP1.h>

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

std::vector<int> getVertsFromPoint(int point, int numVerts, DM &distributed_mesh);

TEST_CASE("Test mapping from reference tet to physical coords and back","[element/tetrahedra]") {

  VectorXd rn,sn,tn;
  std::tie(rn,sn,tn) = Tetrahedra<TetP1>::QuadraturePoints(3);
  int num_pts = rn.size();
  // A distorted element
  Matrix<double,4,3> coord;

  // tet with standard base, but 4th vertex is at (0.5,0.5,0.5) instead of (0,0,1)
  coord <<
    0,0,0,
    0,1,0,    
    1,0,0,    
    0.5,0.5,0.5;
    
  Eigen::VectorXd pts_x,pts_y,pts_z;
  std::tie(pts_x,pts_y,pts_z) = TetP1::buildNodalPoints(rn,sn,tn,coord);

  VectorXd ref_pts_x(num_pts);
  VectorXd ref_pts_y(num_pts);
  VectorXd ref_pts_z(num_pts);

  auto x1 = coord(0,0);
  auto x2 = coord(1,0);
  auto x3 = coord(2,0);
  auto x4 = coord(3,0);
            
  auto y1 = coord(0,1);
  auto y2 = coord(1,1);
  auto y3 = coord(2,1);
  auto y4 = coord(3,1);
            
  auto z1 = coord(0,2);
  auto z2 = coord(1,2);
  auto z3 = coord(2,2);
  auto z4 = coord(3,2);
  
  for (auto i = 0; i < num_pts; i++) {
        
    double r = rn[i];
    double s = sn[i];
    double t = tn[i];
    
    ref_pts_x[i] = r*x3 + s*x2 + t*x4 - x1*(r + s + t - 1);
    ref_pts_y[i] = r*y3 + s*y2 + t*y4 - y1*(r + s + t - 1);
    ref_pts_z[i] = r*z3 + s*z2 + t*z4 - z1*(r + s + t - 1);
  }
  
  REQUIRE((pts_x-ref_pts_x).array().abs().maxCoeff() < 1e-4);
  REQUIRE((pts_y-ref_pts_y).array().abs().maxCoeff() < 1e-4);
  REQUIRE((pts_z-ref_pts_z).array().abs().maxCoeff() < 1e-4);

  
  VectorXd ref_pts_r(num_pts);
  VectorXd ref_pts_s(num_pts);
  VectorXd ref_pts_t(num_pts);
  
  for (auto i = 0; i < num_pts; i++) {
    auto rst = TetP1::inverseCoordinateTransform(pts_x[i],pts_y[i],pts_z[i],coord);
    ref_pts_r[i] = rst[0];
    ref_pts_s[i] = rst[1];
    ref_pts_t[i] = rst[2];
    // printf("(%f,%f,%f) vs (%f,%f,%f)\n",rst[0],rst[1],rst[2],rn[i],sn[i],tn[i]);
  }
  
  // std::cout << "ref_pts_r=\n" << ref_pts_r << "\n";
  
  REQUIRE((rn-ref_pts_r).array().abs().maxCoeff() < 1e-4);
  REQUIRE((sn-ref_pts_s).array().abs().maxCoeff() < 1e-4);
  REQUIRE((tn-ref_pts_t).array().abs().maxCoeff() < 1e-4);
  
}

// TEST_CASE("Testing fluid vs. non-fluid meshes", "[model]") {
  
//   std::cout << "Testing if model crashes without fluid flag.\n";
//   PetscOptionsClear();
//   const char *arg[] = {
//     "salvus_test",
//     "--testing","true",
//     "--duration", "0.08838834764831843", // 120 steps
//     // dt=0.0036084391824351613 is stable for V=1
//     // dt/4 = 0.0009021097956087903
//     "--time_step", "0.0009021097956087903",
//     // "--time_step", "0.0001",
//     "--exodus_file_name", "simple_hexmesh_2x2x2_8elements.vp4.e",
//     "--exodus_model_file_name", "simple_hexmesh_2x2x2_8elements.vp4.e",
//     "--mesh_type", "newmark",
//     "--element_shape", "tet",
//     "--physics_system", "acoustic",
//     "--polynomial_order", "3",
//     "--dirichlet-boundaries", "x0,x1,y0,y1,z0,z1",
//     "--testIC", "true",
//     "--IC-center-x", "0.0",
//     "--IC-center-y", "0.0",
//     "--IC-center-z", "0.0",
//     "--IC-square-side-L", "2",
//     "--saveMovie","false",
//     "--saveFrameEvery","1",
//     "--output_movie_file_name","/scratch/salvus/output_files/movie.h5",    
//     NULL};
//   char **argv = const_cast<char **> (arg);
//   int argc = sizeof(arg) / sizeof(const char *) - 1;
//   PetscOptionsInsert(&argc, &argv, NULL);
  
//   // Set options for exact tests
//   Options options;
//   options.setOptions();
  
//   // Get mesh.
//   Mesh *mesh = Mesh::factory(options);
//   mesh->read(options);
  
//   // Get model.
//   ExodusModel *model = new ExodusModel(options);
//   model->initializeParallel();

//   // Setup reference element.
//   auto reference_element = Element::factory(options);


//   // Setup the dofs on each mesh point.
//   mesh->setupGlobalDof(reference_element->NumDofVtx(),
//                        reference_element->NumDofEdg(),
//                        reference_element->NumDofFac(),
//                        reference_element->NumDofVol(),
//                        reference_element->NumDim(),
//                        model);

// }


// TEST_CASE("Test tetrahedra closure mapping","[element/tetrahedra]") {

//   int order = 3;
  
//   PetscOptionsClear();
//   const char *arg[] = {
//     "salvus_test",
//     "--testing","true",
//     "--exodus_file_name", "simple_tetmesh_2x2x2_30elements.vp4.e",
//     "--exodus_model_file_name", "simple_tetmesh_2x2x2_30elements.vp4.e",
//     "--mesh_type", "newmark",
//     "--element_shape", "tet",
//     "--physics_system", "acoustic",
//     "--polynomial_order", "3",NULL};

//   char **argv = const_cast<char **> (arg);
//   int argc = sizeof(arg) / sizeof(const char *) - 1;
//   PetscOptionsInsert(&argc, &argv, NULL);

//   // Set options for exact tests
//   Options options;
//   options.setOptions();

//   // Setup reference element.
//   auto reference_element = Element::factory(options);
//   auto ref_tet = std::dynamic_pointer_cast<AcousticTet> (reference_element);
//   // Get mesh.
//   Mesh *mesh = Mesh::factory(options);

//     int mesh_load_option = 1;
//   MatrixXd vertices;
//   MatrixXi cells;
//   if(mesh_load_option == 0) {
//     mesh->read(options);  
//   }
//   else if(mesh_load_option == 1) {
//     vertices.resize(3,5);
//     cells.resize(4,2);
//     MatrixXd vertices_t(5,3);
//     vertices_t <<
//       0, 0, 0,
//       0, 1, 0,
//       1, 0, 0,
//       0, 0, 1,
//       0, 0, -1;
//     // -1, 0, 0;
//     vertices = vertices_t.transpose();
//     MatrixXi cells_t(2,4);
        
//     // (-3)
//     // cells_t <<
//     //   0, 1, 2, 3,
//     //   0, 2, 1, 4;
    
//     // (-2)
//     // cells_t <<
//     // 0, 1, 2, 3,
//     // 1, 3, 0, 4;
    
//     // (-1)
//     // cells_t <<
//     // 0, 1, 2, 3,
//     // 3, 0, 1, 4;

//     // bottom to bottom (-3)
//     cells_t <<
//       0, 1, 2, 3,
//       0, 2, 1, 4;
    
//     // bottom to bottom (-2)
//     cells_t <<
//       0, 1, 2, 3,
//       2, 1, 0, 4;
    
//     // bottom to bottom (-1)
//     // cells_t <<
//     // 0, 1, 2, 3,
//     // 1, 0, 2, 4;
    
//     cells = cells_t.transpose();
//     // std::cout << "cells_t=\n" << cells_t << "\n";
    
//     mesh->read(3, 2, 12, 4, cells, vertices);
//     // void Mesh::read(int dim, int numCells, int numVerts, int numVertsPerElem,
//     // Eigen::MatrixXi cells, Eigen::MatrixXd vertex_coords)
//   }
  
//   mesh->setupGlobalDof(ref_tet->NumDofVtx(),
//                        ref_tet->NumDofEdg(),
//                        ref_tet->NumDofFac(),
//                        ref_tet->NumDofVol(),
//                        ref_tet->NumDim());

//   auto Nelem = 1;
  
//   int numPoints;
//   int* points = NULL;
//   DMPlexGetTransitiveClosure(mesh->DistributedMesh(),0,PETSC_TRUE,&numPoints,&points);
//   for(int i=0;i<numPoints;i++) {
//     printf("p[%i]=(%d,%d)\n",i,points[i*2+0],points[i*2+1]);    
//   }

//   auto verts_0 = getVertsFromPoint(1,4,mesh->DistributedMesh());  
//   for(int i=0;i<4;i++) {
//     printf("vert=%d\n",verts_0[i]);
//   }

//   DMPlexGetTransitiveClosure(mesh->DistributedMesh(),1,PETSC_TRUE,&numPoints,&points);
//   for(int i=0;i<numPoints;i++) {
//     printf("p[%i]=(%d,%d)\n",i,points[i*2+0],points[i*2+1]);    
//   }
  
//   std::vector<std::shared_ptr<Element>> elements;  
//   // Get a list of all local elements.
//   for (int i = 0; i < mesh->NumberElementsLocal(); i++) {
//     elements.push_back(std::shared_ptr<Element>(reference_element->clone()));
//   }

//   auto tet0 = std::dynamic_pointer_cast<AcousticTet> (elements[0]);
//   auto tetN = std::dynamic_pointer_cast<AcousticTet> (elements[Nelem]);
  
//   int num_pts = 50;
//   int element_number = 0;
//   for (auto &element_gen : elements) {
//     auto tet = std::dynamic_pointer_cast<AcousticTet> (element_gen);
//     // Give each element a number starting from zero.
//     tet->SetNum(element_number++);
      
//     // Get vertex coordinates from the PETSc DMPLEX.
//     tet->attachVertexCoordinates(mesh->DistributedMesh());          
//   }

//   VectorXd xn,yn,zn;
//   std::tie(xn,yn,zn) = tet0->buildNodalPoints();

//   printf("all points:\n");
//   for(int i=0;i<50;i++) {
//     printf("%d:(%f,%f,%f)\n",i,xn[i],yn[i],zn[i]);
//   }
  
//   for (auto field : mesh->GlobalFields()) {        
//     mesh->registerFieldVectors(field);
//   }
  
//   VectorXd one_at_one_node(num_pts);
//   one_at_one_node.setZero();
//   int selected_node = 16;
//   for(int i=34;i<46;i++) {
//     one_at_one_node[i] = i;
//   }
    
//   mesh->addFieldFromElement("u", 0,
//                             tet0->ClsMap(),
//                             one_at_one_node);

//   one_at_one_node.setZero();
//   one_at_one_node = mesh->getFieldOnElement("u", Nelem,
//                                             tetN->ClsMap());
  
//   for(int i=0;i<num_pts;i++) {
//     printf("%d:%f",i,one_at_one_node[i]);
//     if(one_at_one_node[i] > 0) printf("!!!!!! ----  !!!");
//     printf("\n");
//   }

  
      
// }

