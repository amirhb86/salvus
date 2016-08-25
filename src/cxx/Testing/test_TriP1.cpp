#include "catch.h"
#include <Eigen/Dense>
#include <Utilities/Types.h>
#include <Mesh/Mesh.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Physics/Scalar.h>
#include <Element/ElementAdapter.h>
#include <Problem/Problem.h>
#include <petscviewerhdf5.h>
#include <Element/Simplex/Triangle.h>
#include <Element/Simplex/TriP1.h>

#include <stdexcept>

TEST_CASE("test closure mapping triangle","[tri/closure]") {

  std::string e_file = "tri_eigenfunction.e";

  PetscOptionsClear(NULL);
  const char *arg[] = {
    "salvus_test",
    "--testing", "true",
    "--mesh-file", e_file.c_str(),
    "--model-file", e_file.c_str(),
    "--time-step", "1e-2",
    "--polynomial-order", "3", NULL};
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(NULL, &argc, &argv, NULL);

  std::unique_ptr<Options> options(new Options);
  options->setOptions();

  std::unique_ptr<Problem> problem(Problem::Factory(options));
  std::unique_ptr<ExodusModel> model(new ExodusModel(options));
  std::unique_ptr<Mesh> mesh(Mesh::Factory(options));

  model->read();
  int cells[] =
    {0,1,2,
     1,3,2};
  double vertex_coords[] =
    {0.0,0.0,
     1.0,0.0,
     0.0,1.0,
     1.0,1.0};
     
  mesh->read(2,2,4,3,cells,vertex_coords);

  PetscInt* points=NULL;
  PetscInt numPoints;
  // elemnts are 0,1?
  DMPlexGetTransitiveClosure(mesh->DistributedMesh(),1,PETSC_TRUE,&numPoints,&points);
  // for(int i=0;i<numPoints;i++) {
  //   printf("p[%d]=(%d,o:%d)\n",i,points[2*i],points[2*i+1]);
  // }
  DMPlexRestoreTransitiveClosure(mesh->DistributedMesh(),1,PETSC_TRUE,&numPoints,&points);
     
  /* Setup topology from model and mesh. */
  mesh->setupTopology(model, options);

  /* Setup elements from model and topology. */
  auto elements = problem->initializeElements(mesh, model, options);

  /* Setup global degrees of freedom based on element 0. */
  mesh->setupGlobalDof(elements[0], options);

  std::vector<std::unique_ptr<Element>> test_elements;
  auto fields = problem->initializeGlobalDofs(elements, mesh);
  
  
  RealVec pts_x, pts_y;
  Scalar<Triangle<TriP1>> *e0 = dynamic_cast<Scalar<Triangle<TriP1>>*>(elements[0].get());
  Scalar<Triangle<TriP1>> *e1 = dynamic_cast<Scalar<Triangle<TriP1>>*>(elements[1].get());
  std::tie(pts_x, pts_y) = e0->buildNodalPoints();

  PetscScalar x0 = 0.5, y0 = 0.5, L = 1;
  // RealVec un =
  //   (M_PI / L * (pts_x.array() - (x0 + L / 2))).sin() *
  //   (M_PI/ L * (pts_y.array() - (y0 + L / 2))).sin() *
  //   (M_PI / L * (pts_z.array() - (z0 + L / 2))).sin();
  
  RealVec un =
    (M_PI / L * (pts_x.array() - (x0 + L / 2))).sin();
  
  // for(int i=0;i<un.size();i++) {
  //   un[i] = pts_y[i];
  // }
  
  RealMat grad_xyz0 = e0->computeGradient(un);
  RealMat grad_xyz1= e1->computeGradient(un);
  
  // problem->insertElementalFieldIntoMesh("u", 0, elements[0]->ClsMap(), un,
  //                                       mesh->DistributedMesh(), mesh->MeshSection(),
  //                                       fields);

  un.setZero(pts_x.size());
  for(int i=0;i<12;i++) {
    un[i] = i;
  }
  problem->addFieldOnElement("u",0,elements[0]->ClsMap(), un,
                             mesh->DistributedMesh(), mesh->MeshSection(),
                             fields);

  // un.setZero(pts_x.size());
  // un[5] = 1.0;
  // problem->addFieldOnElement("u",1,elements[1]->ClsMap(), un,
  //                            mesh->DistributedMesh(), mesh->MeshSection(),
  //                            fields);
  
  RealVec u_1 = problem->getFieldOnElement("u", 1, elements[1]->ClsMap(),
                                         mesh->DistributedMesh(), mesh->MeshSection(), fields);
  
  // for(int i=0;i<16;i++) {
  //   if(u_1[i] > 0.0) {
  //     printf("u_1[%d]=%f\n",i,u_1[i]);
  //   }
  // }
  
  std::vector<int> ids = {11,7,8,9};
  std::vector<double> ans = {11.0,6.0,5.0,10.0};
  for(int i=0;i<ids.size();i++) {
    REQUIRE(u_1[ids[i]] == ans[i]);
  }
  
  
  // for(int i=0;i<64;i++) {
  // printf("grad_x0[%d]=%f\n",i,grad_xyz0(i,0));
  // }
  // for(int i=0;i<64;i++) {
  // printf("grad_x1[%d]=%f\n",i,grad_xyz1(i,0));
  // }
  
  
}
