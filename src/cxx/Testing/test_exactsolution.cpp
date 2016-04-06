#include "catch.h"
#include <Eigen/Dense>
#include <petsc.h>
#include <Element/Element.h>
#include <Element/HyperCube/Quad/Acoustic.h>
#include <Element/Simplex/Triangle/AcousticTri.h>

int initialize_exact(Mesh *mesh,
                     ExodusModel *model,
                     Element *reference_element,
                     Options options, std::vector<Element*>* elements) {
    
    PetscFunctionBegin;

    // Attach elements to mesh.
    mesh->setupGlobalDof(reference_element->NumDofVtx(),
                         reference_element->NumDofEdg(),
                         reference_element->NumDofFac(),
                          0 /* zero dofvolume */,
                         reference_element->NumDim());

    // Setup boundary conditions from options.
    mesh->setupBoundaries(options);

    // Register all global field required for time stepping.
    for (auto field : mesh->GlobalFields()) {        
        mesh->registerFieldVectors(field);
    }

    // Get a list of all local elements.
    for (int i = 0; i < mesh->NumberElementsLocal(); i++) {
        elements->push_back(reference_element->clone());
    }    
    
    // Set up elements.
    int element_number = 0;
    for (auto &element : *elements) {

        // Give each element a number starting from zero.
        element->SetNum(element_number++);

        // Get vertex coordinates from the PETSc DMPLEX.
        element->attachVertexCoordinates(mesh->DistributedMesh());

        // Add material parameters (velocity, Cij, etc...).
        element->attachMaterialProperties(model);

        // Set boundary conditions.
        element->setBoundaryConditions(mesh);

        // Assemble the (elemental) mass matrix.
        element->assembleElementMassMatrix(mesh);

        // Prepare stiffness terms
        element->prepareStiffness();

        // setup tests
        element->setupTest(mesh, options);
        
    }

    // Scatter the mass matrix to the global dofs.
    mesh->assembleLocalFieldToGlobal("m");

    // set initial condition on global vectors
    for (auto &field : reference_element->PullElementalFields()) {
        mesh->setLocalFieldToGlobal(field);
    }

    PetscFunctionReturn(0);
    
}

double solve_vs_exact(Options options, Mesh* mesh, std::vector<Element*> &elements) {
    PetscFunctionBegin;
    // Setup values.
    int it          = 0;
    double time     = 0.0;
    double timeStep = options.TimeStep();
    double duration = options.Duration();
    
    if(options.SaveMovie()) mesh->setUpMovie(options.OutputMovieFile());

    // Over-allocate matrices to avoid re-allocations.
    int max_dims = 3;
    int int_pnts = elements[0]->NumIntPnt();
    Eigen::MatrixXd f(int_pnts, max_dims);
    Eigen::MatrixXd u(int_pnts, max_dims);
    Eigen::MatrixXd ku(int_pnts, max_dims);
    Eigen::MatrixXd fMinusKu(int_pnts, max_dims);
    Eigen::MatrixXd fMinusKu_boundary(int_pnts, max_dims);

    double max_error = 0.0;    

    while (time < duration) {

        // Collect all global fields to the local partitions.
        for (auto &field : elements[0]->PullElementalFields()) {
            mesh->checkOutField(field);
        }

        // Zero all fields to which we will sum.
        for (auto &field : elements[0]->PushElementalFields()) {
            mesh->zeroFields(field);
        }

        for (auto &element : elements) {

            // Get fields on element, store in successive rows. (e.g.,
            // "u" for acoustic, "ux, uz" for elastic)
            int fitr = 0;
            for (auto &field : element->PullElementalFields()) {
                u.col(fitr) = mesh->getFieldOnElement(field, element->Num(),
                                                      element->ClsMap());
                fitr++;
            }
            
            double element_error = element->checkTest(mesh, options, u.block(0,0,int_pnts,fitr), time);
            
            if(element_error > max_error) { max_error = element_error; }

            // Compute stiffness, only passing those rows which are occupied.
            ku.block(0,0,int_pnts,fitr) = element->computeStiffnessTerm(
                                                                        u.block(0,0,int_pnts,fitr));

            // Compute acceleration.
            fMinusKu.block(0,0,int_pnts,fitr) = -1 * ku.block(0,0,int_pnts,fitr).array();

            
            // Sum fields into local partition.
            fitr = 0;
            for (auto &field : element->PushElementalFields()) {
                mesh->addFieldFromElement(field, element->Num(),
                                          element->ClsMap(),
                                           fMinusKu.col(fitr));
                fitr++;
            }
        }

        // boundary condition happens after full Ku
        for( auto &element : elements) {
            // we can now apply boundary conditions (after fields are set)
            if(element->BndElm()) {
                // apply boundary condition
                element->applyBoundaryConditions(mesh,
                                                 options,
                                                 "a");
            }
        }


        // Sum fields into global partitions.
        for (auto &field : elements[0]->PushElementalFields()) {
            mesh->checkInFieldBegin(field);
            mesh->checkInFieldEnd(field);
        }

        // Take a time step.
        mesh->applyInverseMassMatrix();
        mesh->advanceField(timeStep);
        
        if(options.SaveMovie() && (it%options.SaveFrameEvery()==0 || it == 0) ) {
            // GlobalFields[0] == "u" for acoustic and == "ux" for elastic
            mesh->saveFrame(mesh->GlobalFields()[0], it);
            // mesh->setLocalFieldToGlobal("u_exact");
            // mesh->saveFrame("u_exact", it);
            if(max_error > 5) {
                std::cerr << "ERROR: Solution blowing up!\n";
                exit(1);
            }
            PRINT_ROOT() << "TIME: " << time;
        }
        
        it++;
        time += timeStep;

    }

    PRINT_ROOT() << "Max Error @ T=end: " << max_error << std::endl;

    if(options.SaveMovie()) mesh->finalizeMovie();
    return max_error;
}

TEST_CASE("Testing acoustic exact solutions for triangles", "[exact/triangles]") {

    // Set options for exact tests
    Options options;
    // Triangles
    options.__SetPolynomialOrder(3);
    options.__SetDuration(0.7071067811865475); // 1 full cycle with v=4
    options.__SetTimeStep(0.003);
    options.__SetMeshType("newmark");
    options.__SetExodusMeshFile("simple_trimesh_2x2.e");
    options.__SetExodusModelFile("simple_trimesh_2x2.e");
    options.__SetElementShape("triangle");
    options.__SetPhysicsSystem("acoustic");
    options.__SetDirichletBoundaryNames({"dirichlet"});
    options.__SetCenter_x(0.0);
    options.__SetCenter_z(0.0);
    options.__SetSquareSide_L(2.0);
    options.__SetSaveMovie(PETSC_FALSE);
    
    
    // Get mesh.
    Mesh *mesh = Mesh::factory(options);
    mesh->read(options);
    
    // Get model.
    ExodusModel *model = new ExodusModel(options);
    model->initializeParallel();

    // Setup reference element.
    Element *reference_element = Element::factory(options);
    
    std::vector<Element*> elements;

    initialize_exact(mesh, model, reference_element, options, &elements);
    
    double error = solve_vs_exact(options, mesh, elements);
    REQUIRE(error < 3e-2);
    
}

TEST_CASE("Testing acoustic exact solutions for quadrilaterals", "[exact/quads]") {

    // Set options for exact tests
    Options options;
    // Triangles
    options.__SetPolynomialOrder(3);
    options.__SetDuration(0.7071067811865475); // 1 full cycle with v=4
    options.__SetTimeStep(0.003);
    options.__SetMeshType("newmark");
    options.__SetExodusMeshFile("simple_quadmesh_2x2.e");
    options.__SetExodusModelFile("simple_quadmesh_2x2.e");
    options.__SetElementShape("quad");
    options.__SetPhysicsSystem("acoustic");
    options.__SetDirichletBoundaryNames({"x0","x1","y0","y1"});
    options.__SetCenter_x(0.0);
    options.__SetCenter_z(0.0);
    options.__SetSquareSide_L(2.0);
    options.__SetSaveMovie(PETSC_FALSE);
    
    // Get mesh.
    Mesh *mesh = Mesh::factory(options);
    mesh->read(options);
    
    // Get model.
    ExodusModel *model = new ExodusModel(options);
    model->initializeParallel();

    // Setup reference element.
    Element *reference_element = Element::factory(options);
    
    std::vector<Element*> elements;

    initialize_exact(mesh, model, reference_element, options, &elements);
    
    double error = solve_vs_exact(options, mesh, elements);
    REQUIRE(error < 3e-2);
    
}
