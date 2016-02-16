//
// Created by Michael Afanasiev on 2016-01-27.
//

#include "Problem.h"


Problem *Problem::factory(std::string solver_type) {
    try {

        if (solver_type == "time_domain") {
            return new TimeDomain;
        } else {
            throw std::runtime_error("Runtime Error: Problem type " + solver_type + " not supported.");
        }
    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI::COMM_WORLD.Abort(-1);
        return nullptr;
    }
}

void TimeDomain::initialize(Mesh *mesh, ExodusModel *model, Quad *quad, Options options) {

    // Perform dynamic casts to ensure types are appropriate.
    mMesh = dynamic_cast<ScalarNewmark*> (mesh);
    mReferenceQuad = dynamic_cast<Quad*> (quad);

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
        element->attachVertexCoordinates(mesh->DistributedMesh());
        element->interpolateMaterialProperties(model);
        element->readOperators();
        element->assembleMassMatrix();
        element->scatterMassMatrix(mesh);
    }

    // Scatter the mass matrix to the global dof.
    mesh->checkInFieldBegin("mass_matrix");
    mesh->checkInFieldEnd("mass_matrix");

    // Set up options.
    mMesh->setUpMovie(options.OutputMovieFile());
    mSimulationDuration = options.Duration();
    mTimeStep = options.TimeStep();

}

void TimeDomain::solve() {

    double time = 0;
    while (time < mSimulationDuration) {

        // Pull down the displacement (pressure) from the global dof.
        mMesh->checkOutField("displacement");

        // Compute element-wise terms.
        for (auto &element: mElements) {

            // Check out the necessary fields on each element (i.e. displacement) from the mesh.
            element->checkOutField(mMesh);

            // Set the time on the element (useful for things like source firing).
            element->SetTime(time);

            // Compute the stiffness term (apply K, if you will).
            element->computeStiffnessTerm();

            // Fire any sources and handle the different cases.
            element->computeSourceTerm();

            // Compute any surface integral terms (i.e. if we're in a coupling layer).
            element->computeSurfaceTerm();

            // Sum the element fields (i.e. forcing) back into the locally-owned section of the mesh.
            element->checkInField(mMesh);

        }

        // Sum fields back into global dof.
        mMesh->checkInFieldBegin("force");
        mMesh->checkInFieldEnd("force");

        // Timestep specific things.
        mMesh->applyInverseMassMatrix();
        mMesh->advanceField();

        // Save to file.
        mMesh->saveFrame();

        std::cout << time << std::endl;
        time += mTimeStep;

    }

    mMesh->finalizeMovie();

}


