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
//        Eigen::VectorXd pts_x,pts_z;
//        std::tie(pts_x,pts_z) = element->buildNodalPoints(mesh);
//        element->setInitialCondition(mesh,pts_x,pts_z);
    }

    // Assemble the mass matrix to the global dof.
    mesh->assembleLocalFieldToGlobal("mass_matrix");
//    mesh->setLocalFieldToGlobal("nodes_x");
//    mesh->setLocalFieldToGlobal("nodes_z");

    // put initial condition on global dofs
    mesh->setLocalFieldToGlobal("displacement");

    // Set up options.
    mMesh->setUpMovie(options.OutputMovieFile());
    mSimulationDuration = options.Duration();
    mTimeStep = options.TimeStep();

}

void TimeDomainScalar2d::solve() {

    Eigen::MatrixXd F;
    Eigen::MatrixXd u;
    Eigen::MatrixXd Ku;
    Eigen::MatrixXd FminusKu;

    double time = 0;
    while (time < mSimulationDuration) {

        // Pull down the displacement (pressure) from the global dof.
        mMesh->checkOutField("displacement");
        mMesh->zeroFields("force");

        // Compute element-wise terms.
        for (auto &element: mElements) {

            // Set the time on the element (useful for things like source firing).
            element->SetTime(time);

            // Check out the necessary fields on each element (i.e. displacement) from the mesh.
            u = element->checkOutFieldElement(mMesh, "displacement");

            // Compute the stiffness term (apply K, if you will).
            Ku = element->computeStiffnessTerm(u);

            // Fire any sources and handle the different cases.
            F = element->computeSourceTerm();

            // Compute the force balance.
            FminusKu = F - Ku;

            // Compute any surface integral terms (i.e. if we're in a coupling layer).
            element->computeSurfaceTerm();

            // Sum the element fields (i.e. forcing) back into the locally-owned section of the mesh.
            element->checkInFieldElement(mMesh,FminusKu,"force");

        }


        // Sum fields back into global dof.
        mMesh->checkInFieldBegin("force");
        mMesh->checkInFieldEnd("force");

        // Time step specific things.
        mMesh->applyInverseMassMatrix();
        mMesh->advanceField();

        // Save to file.
        mMesh->saveFrame("displacement");

        std::cout << time << std::endl;
        time += mTimeStep;

    }

    mMesh->finalizeMovie();

}