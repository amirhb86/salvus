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
    mMesh->setupBoundaries(options);

    // Initialize all fields relevant to the chosen time-stepping scheme.
    mMesh->registerFields();

    // Clone a list of elements to do computations on.
    for (auto i = 0; i < mMesh->NumberElementsLocal(); i++) {
        mElements.push_back(mReferenceQuad->clone());
    }

    // Set up the properties on the individual elements (material parameters, sources, etc.)
    int element_number = 0;
    for (auto &element: mElements) {
        element->SetLocalElementNumber(element_number++);
        element->attachVertexCoordinates(mesh->DistributedMesh());
        element->interpolateMaterialProperties(model);
        element->setBoundaryConditions(mesh);
        element->readGradientOperator();
        element->assembleElementMassMatrix(mesh);
        element->attachSource(sources);
//        Eigen::VectorXd pts_x,pts_z;
//        std::tie(pts_x,pts_z) = element->buildNodalPoints(mesh);
//        element->setInitialCondition(mesh,pts_x,pts_z);
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
    double max_error = 0.0;
    while (time < mSimulationDuration) {

        // Pull down the displacement (pressure) from the global dof.
        mMesh->checkOutField("displacement");
        mMesh->zeroFields("force");

        max_error = 0.0;
        
        // Compute element-wise terms.
        for (auto &element: mElements) {

            // Check out the necessary fields on each element (i.e. displacement) from the mesh.
            u = mMesh->getFieldOnElement("displacement",
                                         element->Number(),
                                         element->ElementClosure());

            // exact solution
//            auto x_e = element->checkOutFieldElement(mMesh, "nodes_x");
//            auto z_e = element->checkOutFieldElement(mMesh, "nodes_z");
//            auto un_exact = element->exactSolution(x_e,z_e,time);
//            element->setFieldElement(mMesh,un_exact,"displacement_exact");
            
//            auto element_error = (un_exact - u).array().abs().maxCoeff();
//            if(element_error > max_error) { max_error = element_error; }

            // Compute the stiffness term (apply K, if you will).
            Ku = element->computeStiffnessTerm(u);

            // Fire any sources and handle the different cases.
            F = element->computeSourceTerm(time);

            // Compute the force balance.
            FminusKu = F - Ku;

            // Sum the element fields back into the locally-owned section of the mesh.
            mMesh->addFieldFromElement("force",
                                       element->Number(),
                                       element->ElementClosure(),
                                       FminusKu);

        }

        // Sum fields back into global dof.
        mMesh->checkInFieldBegin("force");
        mMesh->checkInFieldEnd("force");

        // Once all element-wise terms are calculated, go back and compute boundary terms.
        for (auto &element: mElements) {

            for (auto &boundary: element->Boundaries()) {

                auto b_val = mMesh->getFieldOnFace("force",
                                                   boundary.second,
                                                   element->GetFaceClosureMapping());
                element->computeSurfaceTerm();

            }

        }

        exit(-1);

        applyDirichletBoundary(mMesh,"force",0.0);        

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
        
        if(it%options.SaveFrameEvery()==0) {
            if(max_error > 5) {
                std::cerr << "ERROR: Solution blowing up!\n";
                exit(1);
            }
            PRINT_ROOT() << "Time(" << it << "): " << time;
            std::cout << "Max Error: " << max_error << std::endl;
        }
        
        time += mTimeStep;
        it++;
    }
    PRINT_ROOT() << "Max Error @ T=end: " << max_error << std::endl;

    if(options.SaveMovie()) { mMesh->finalizeMovie(); }

}


void TimeDomainScalar2d::applyDirichletBoundary(Mesh *mesh,std::string fieldname, double value) {

//    auto& boundaryIds = mesh->BoundaryIds();
//    auto dirichlet_id_it = boundaryIds.find("dirichlet");
//    // found dirichlet boundary
//    if(dirichlet_id_it != boundaryIds.end()) {
//        int dirichlet_id = (*dirichlet_id_it).second;
//        auto& dirichlet_elems = mesh->BoundaryElementFaces()[dirichlet_id];
//        for(auto& elem_key : dirichlet_elems) {
//            for(auto& faceid : elem_key.second) {
//                auto field = mesh->getFieldOnFace(fieldname,faceid,
//                                                   mElements[elem_key.first]->GetFaceClosureMapping());
//                // apply dirichlet condition
//                field = 0*field.array() + value;
//                mesh->setFieldFromFace(fieldname,faceid,mElements[elem_key.first]->GetFaceClosureMapping(),field);
//            }
//        }
//    }
}
