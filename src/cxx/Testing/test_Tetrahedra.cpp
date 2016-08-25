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
#include <Element/Simplex/Tetrahedra.h>
#include <Element/Simplex/TetP1.h>

#include <stdexcept>

TEST_CASE("test closure mapping tetrahedra","[tet/closure]") {

  std::string e_file = "tet_eigenfunction.e";

  std::vector<std::vector<int>> triangle2 = {{0,2,1,4}, // o:-3
                                             {2,1,0,4}, // o:-2
                                             {1,0,2,4} }; // o:-1

  for(int n=0;n<triangle2.size();n++) {
  
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
      {0,1,2,3,
       triangle2[n][0],triangle2[n][1],triangle2[n][2],triangle2[n][3]
      };
    double vertex_coords[] =
      {0.0,0.0,0.0,
       0.0,1.0,0.0,
       1.0,0.0,0.0,
       0.0,0.0,1.0,
       0.0,0.0,-1.0
      };
     
    mesh->read(3,2,5,4,cells,vertex_coords);

    PetscInt* points=NULL;
    PetscInt numPoints;

    // elemnts are 0,1?
    DMPlexGetTransitiveClosure(mesh->DistributedMesh(),0,PETSC_TRUE,&numPoints,&points);
    // for(int i=0;i<numPoints;i++) {
    //   printf("p[%d]=(%d,o:%d)\n",i,points[2*i],points[2*i+1]);
    // }
    DMPlexRestoreTransitiveClosure(mesh->DistributedMesh(),0,PETSC_TRUE,&numPoints,&points);
    // printf("\n");
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
  
    RealVec pts_x, pts_y, pts_z;
    Scalar<Tetrahedra<TetP1>> *e0 = dynamic_cast<Scalar<Tetrahedra<TetP1>>*>(elements[0].get());
    Scalar<Tetrahedra<TetP1>> *e1 = dynamic_cast<Scalar<Tetrahedra<TetP1>>*>(elements[1].get());
    std::tie(pts_x, pts_y, pts_z) = e0->buildNodalPoints();

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
    for(int i=0;i<50;i++) {
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
  
    // for(int i=0;i<50;i++) {
    //   if(u_1[i] > 0.0) {
    //     printf("u_1[%d]=%f\n",i,u_1[i]);
    //   }
    // }
    
    std::vector<int> ids = {10,11,12,13,14,15,
                            34,35,36,37,38,39,
                            46,47,48};
    std::vector<double> ans = {10,12,11,13,15,14,
                               39,38,37,36,35,34,
                               46,48,47};
    std::vector<std::vector<double>> answers = {{10,12,11,13,15,14,
                                                 39,38,37,36,35,34,
                                                 46,48,47},
                                                {12,11,10,15,14,13,
                                                 37,36,35,34,39,38,
                                                 48,47,46},
                                                {11,10,12,14,13,15,
                                                 35,34,39,38,37,36,
                                                 47,46,48}};
    for(int i=0;i<ids.size();i++) {
      REQUIRE(u_1[ids[i]] == answers[n][i]);
    }
  
  }
  // for(int i=0;i<64;i++) {
  // printf("grad_x0[%d]=%f\n",i,grad_xyz0(i,0));
  // }
  // for(int i=0;i<64;i++) {
  // printf("grad_x1[%d]=%f\n",i,grad_xyz1(i,0));
  // }
  
  
}
