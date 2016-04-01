#include "NewmarkTesting.h"

void NewmarkTesting::initialize(Mesh *mesh,
                                ExodusModel *model,
                                Element2D *elem,
                                Options options) {

    // Save references to mesh and element base.
    mMesh = mesh;
    mReferenceElem = elem;

    // Attach elements to mesh.
    mMesh->setupGlobalDof(mReferenceElem->NumberDofVertex(),
                          mReferenceElem->NumberDofEdge(),
                          mReferenceElem->NumberDofFace(),
                          0 /* zero dofvolume */,
                          mReferenceElem->NumberDimensions());

    // Setup boundary conditions from options.
    mMesh->setupBoundaries(options);

    // Register all global field required for time stepping.
    for (auto field : mMesh->GlobalFields()) {        
        mMesh->registerFieldVectors(field);
    }

    // Get a list of all local elements.
    for (int i = 0; i < mMesh->NumberElementsLocal(); i++) {
        mElements.push_back(mReferenceElem->clone());
    }

    // Get a list of all sources.
    auto sources = Source::factory(options);

    // Set up elements.
    int element_number = 0;
    for (auto &element : mElements) {

        // Give each element a number starting from zero.
        element->SetLocalElementNumber(element_number++);

        // Get vertex coordinates from the PETSc DMPLEX.
        element->attachVertexCoordinates(mMesh->DistributedMesh());

        // Add material parameters (velocity, Cij, etc...).
        element->interpolateMaterialProperties(model);

        // Set boundary conditions.
        element->setBoundaryConditions(mMesh);

        // Assemble the (elemental) mass matrix.
        element->assembleElementMassMatrix(mMesh);

        // Attach any external source terms.
        element->attachSource(sources);
        
        // Prepare stiffness terms
        element->prepareStiffness();

        // setup tests
        element->setupTest(mMesh, options);
        
    }
    
    // Scatter the mass matrix to the global dofs.
    mMesh->assembleLocalFieldToGlobal("m");

    // set initial condition on global vectors
    for (auto &field : mReferenceElem->PullElementalFields()) {
        mMesh->setLocalFieldToGlobal(field);
    }
    // set nodal locations on global dof (for 2D x,z & 3D x,y,z)
    mMesh->setLocalFieldToGlobal("x");
    if(mMesh->NumberDimensions() == 3) {
        mMesh->setLocalFieldToGlobal("y");
    }
    mMesh->setLocalFieldToGlobal("z");
    
}

void NewmarkTesting::solve(Options options) {

    // Setup values.
    int it          = 0;
    double time     = 0.0;
    double timeStep = options.TimeStep();
    double duration = options.Duration();
    mMesh->setUpMovie(options.OutputMovieFile());

    // Over-allocate matrices to avoid re-allocations.
    int max_dims = 3;
    int int_pnts = mReferenceElem->NumberIntegrationPoints();
    Eigen::MatrixXd f(int_pnts, max_dims);
    Eigen::MatrixXd u(int_pnts, max_dims);
    Eigen::MatrixXd ku(int_pnts, max_dims);
    Eigen::MatrixXd fMinusKu(int_pnts, max_dims);
    Eigen::MatrixXd fMinusKu_boundary(int_pnts, max_dims);

    double max_error = 0.0;
    
    while (time < duration) {

        // Collect all global fields to the local partitions.
        for (auto &field : mReferenceElem->PullElementalFields()) {
            mMesh->checkOutField(field);
        }

        // Zero all fields to which we will sum.
        for (auto &field : mReferenceElem->PushElementalFields()) {
            mMesh->zeroFields(field);
        }

        for (auto &element : mElements) {

            // Get fields on element, store in successive rows. (e.g.,
            // "u" for acoustic, "ux, uz" for elastic)
            int fitr = 0;
            for (auto &field : mReferenceElem->PullElementalFields()) {
                u.col(fitr) = mMesh->getFieldOnElement(field, element->Number(),
                                                       element->ElementClosure());
                fitr++;
            }
            
            double element_error = element->checkTest(mMesh, options, u.block(0,0,int_pnts,fitr), time);
            
            if(element_error > max_error) { max_error = element_error; }
                       
            // Compute stiffness, only passing those rows which are occupied.
            ku.block(0,0,int_pnts,fitr) = element->computeStiffnessTerm(
                    u.block(0,0,int_pnts,fitr));

            // Compute source term.
            f.block(0,0,int_pnts,fitr) = element->computeSourceTerm(time);

            // Compute acceleration.
            fMinusKu.block(0,0,int_pnts,fitr) = f.block(0,0,int_pnts,fitr).array() -
                    ku.block(0,0,int_pnts,fitr).array();

            
            // Sum fields into local partition.
            fitr = 0;
            for (auto &field : mReferenceElem->PushElementalFields()) {
                mMesh->addFieldFromElement(field, element->Number(),
                                           element->ElementClosure(),
                                           fMinusKu.col(fitr));
                fitr++;
            }

        }

        // boundary condition happens after full Ku
        for( auto &element : mElements) {
            // we can now apply boundary conditions (after fields are set)
            if(element->OnBoundary()) {
                // apply boundary condition
                element->applyBoundaryConditions(mMesh,
                                                 options,
                                                 "a");
            }
        }

        // Sum fields into global partitions.
        for (auto &field : mReferenceElem->PushElementalFields()) {
            mMesh->checkInFieldBegin(field);
            mMesh->checkInFieldEnd(field);
        }

        // Take a time step.
        mMesh->applyInverseMassMatrix();
        mMesh->advanceField(timeStep);
        
        if(options.SaveMovie() && (it%options.SaveFrameEvery()==0 || it == 0) ) {
            // GlobalFields[0] == "u" for acoustic and == "ux" for elastic
            mMesh->saveFrame(mMesh->GlobalFields()[0], it);
            // mMesh->setLocalFieldToGlobal("u_exact");
            // mMesh->saveFrame("u_exact", it);
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

    if(options.SaveMovie()) mMesh->finalizeMovie();

}
