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

void TimeDomain::initialize(Mesh *mesh, ExodusModel *model, Quad *quad) {

    // Perform dynamic casts to ensure types are appropriate.
    mMesh = dynamic_cast<ScalarNewmark*> (mesh);
    mReferenceQuad = dynamic_cast<Quad*> (quad);

    // Attach element to mesh.
    mMesh->setupGlobalDof(quad->NumberDofVertex(), quad->NumberDofEdge(),
                          quad->NumberDofFace(), quad->NumberDofVolume(),
                          quad->NumberDimensions());

    // Initialize all fields relevant to the chosen time-stepping scheme.
    mMesh->registerFields();

    // Clone a list of elements to do computations on.
    std::vector<Quad *> elements;
    for (auto i = 0; i < mMesh->NumberElementsLocal(); i++) {
        elements.push_back(mReferenceQuad->clone());
    }

    // Set up the properties on the individual elements (material parameters, sources, etc.)
    int element_number = 0;
    for (auto &element: elements) {
        element->SetLocalElementNumber(element_number);
        element->attachVertexCoordinates(mesh->DistributedMesh());
//        element->interpolateMaterialProperties(model);
    }







}

void TimeDomain::solve() { }
