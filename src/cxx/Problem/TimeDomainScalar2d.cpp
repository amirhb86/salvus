//
// Created by Michael Afanasiev on 2016-02-23.
//

#include "TimeDomainScalar2d.h"

void TimeDomainScalar2d::initialize(Mesh *mesh, ExodusModel *model, Quad *quad, Options options) {

    std::vector<Source *> sources = Source::factory(options);

    // Perform dynamic casts to ensure types are appropriate.
    mMesh = dynamic_cast<ScalarNewmark2D *> (mesh);
    mReferenceQuad = dynamic_cast<Acoustic *> (quad);

    // Attach element to mesh.
    mMesh->setupGlobalDof(mReferenceQuad->NumberDofVertex(), mReferenceQuad->NumberDofEdge(),
                          mReferenceQuad->NumberDofFace(), mReferenceQuad->NumberDofVolume(),
                          mReferenceQuad->NumberDimensions());

    // Initialize all fields relevant to the chosen time-stepping scheme.
    mMesh->registerFields();

    // Clone a list of elements to do computations on.
    for (auto i = 0; i < mMesh->NumberElementsLocal(); i++) {
        mElements.push_back(mReferenceQuad->clone());
    }

    // Set up the properties on the individual elements (material parameters, sources, etc.)
    int element_number = 0;
    for (auto &element: mElements) {
        element->SetLocalElementNumber(element_number);
        element_number++;

        element->attachVertexCoordinates(mesh->DistributedMesh());
        element->interpolateMaterialProperties(model);
        element->readGradientOperator();
        element->assembleElementMassMatrix(mesh);
        element->attachSource(sources);
        Eigen::VectorXd pts_x,pts_z;
        std::tie(pts_x,pts_z) = element->buildNodalPoints(mesh);
        element->setInitialCondition(mesh,pts_x,pts_z);
    }    
    
    // Assemble the mass matrix to the global dof.
    mesh->assembleLocalFieldToGlobal("mass_matrix");
    mesh->setLocalFieldToGlobal("nodes_x");
    mesh->setLocalFieldToGlobal("nodes_z");
    
    // put initial condition on global dofs
    mesh->setLocalFieldToGlobal("displacement");

    mSimulationDuration = options.Duration();
    mTimeStep = options.TimeStep();

    if(options.SaveMovie()) {
        // Set up options.
        mMesh->setUpMovie(options.OutputMovieFile());        
    }

}

void TimeDomainScalar2d::solve(Options options) {

    Eigen::MatrixXd F;
    Eigen::MatrixXd u;
    Eigen::MatrixXd Ku;
    Eigen::MatrixXd FminusKu;

    double time = 0;
    int it=0;
    std::cout << "mSimulationDuration=" << mSimulationDuration << "\n";
    while (time < mSimulationDuration) {

        // Pull down the displacement (pressure) from the global dof.
        mMesh->checkOutField("displacement");
        mMesh->zeroFields("force");

        double max_error = 0.0;
        
        // Compute element-wise terms.
        for (auto &element: mElements) {

            // Set the time on the element (useful for things like source firing).
            element->SetTime(time);

            // Check out the necessary fields on each element (i.e. displacement) from the mesh.
            u = element->checkOutFieldElement(mMesh, "displacement");

            // exact solution
            auto x_e = element->checkOutFieldElement(mMesh, "nodes_x");
            auto z_e = element->checkOutFieldElement(mMesh, "nodes_z");
            auto un_exact = element->exactSolution(x_e,z_e,time);
            element->setFieldElement(mMesh,un_exact,"displacement_exact");
            
            auto element_error = (un_exact - u).array().abs().maxCoeff();
            if(element_error > max_error) { max_error = element_error; }

            // Compute the stiffness term (apply K, if you will).
            Ku = element->computeStiffnessTerm(u);

            // Fire any sources and handle the different cases.
            F = element->computeSourceTerm();

            // Compute the force balance.
            FminusKu = F - Ku;

            // Compute any surface integral terms (i.e. if we're in a coupling layer).
            // element->computeSurfaceTerm();
            
            // element->computeBoundaryConditions(mMesh,FminusKu);            
            // std::cout << "Max force: " << FminusKu.maxCoeff() << std::endl;
            // Sum the element fields back into the locally-owned section of the mesh.
            element->checkInFieldElement(mMesh,FminusKu,"force");            
                       
        }

        if(max_error > 5) {
            std::cerr << "ERROR: Solution blowing up!\n";
            exit(1);
        }
        
        auto& boundaryIds = mMesh->BoundaryIds();
        auto dirichlet_id_it = boundaryIds.find("dirichlet");
        int myrank;
        MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
        // found dirichlet boundary
        if(dirichlet_id_it != boundaryIds.end()) {
            int dirichlet_id = (*dirichlet_id_it).second;
            // printf("dirichlet id=%d\n",dirichlet_id);
            
            auto& dirichlet_elems = mMesh->BoundaryElementFaces()[dirichlet_id];
            // printf("Size of dirichlet elems=%d\n",dirichlet_elems.size());
            for(auto& elem_key : dirichlet_elems) {                
                if(myrank == 0) {
                    auto vertices = mElements[elem_key.first]->GetVertexCoordinates();                    

                    // printf("elem(%d)={",elem_key.first);
                    // for(int i=0;i<4;i++) {
                    //     printf("(%f,%f),",vertices.row(0)[i],vertices.row(1)[i]);
                    // }
                    // printf("}\n");
                    
                }
                for(auto& faceid : elem_key.second) {
                    // if(myrank == 0) printf("faceid=%d\n",faceid);
                    // auto nodes_x = mMesh->getFieldOnFace("nodes_x",faceid,
                    //                                      mElements[elem_key.first]->GetFaceClosureMapping());
                    // auto nodes_z = mMesh->getFieldOnFace("nodes_z",faceid,
                    //                                      mElements[elem_key.first]->GetFaceClosureMapping());
                    // for(int xz=0;xz<nodes_x.size();xz++) {
                    // if(myrank == 0) printf("face pt: (%f,%f)\n",nodes_x[xz],nodes_z[xz]);
                    // }
                                        
                    auto field = mMesh->getFieldOnFace("force",faceid,
                                                       mElements[elem_key.first]->GetFaceClosureMapping());
                    // apply 0-dirichlet condition
                    field = 0*field;
                    mMesh->setFieldFromFace("force",faceid,mElements[elem_key.first]->GetFaceClosureMapping(),field);
                    
                    
                }
            }
        }
        
        std::cout << "Max Error: " << max_error << std::endl;
        
        // Sum fields back into global dof.
        mMesh->checkInFieldBegin("force");
        mMesh->checkInFieldEnd("force");

        mMesh->setLocalFieldToGlobal("displacement_exact");
        
        // Time step specific things.
        mMesh->applyInverseMassMatrix();
        mMesh->advanceField(mTimeStep);

        if(options.SaveMovie() && (it%options.SaveFrameEvery()==0 || it == 0)) {
            // Save to file.
            mMesh->saveFrame("displacement",it);
            mMesh->saveFrame("displacement_exact",it);
            // mMesh->saveFrame("acceleration",it);
            // mMesh->saveFrame("velocity",it);
            // mMesh->saveFrame("force",it);
        }        
        
        PRINT_ROOT() << "Time(" << it << "): " << time;
        
        time += mTimeStep;
        it++;

    }

    if(options.SaveMovie()) { mMesh->finalizeMovie(); }

}
