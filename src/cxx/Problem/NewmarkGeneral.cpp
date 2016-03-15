//
// Created by Michael Afanasiev on 2016-03-14.
//

#include "NewmarkGeneral.h"

void NewmarkGeneral::initialize(Mesh *mesh,
                                ExodusModel *model,
                                Quad *quad,
                                Options options) {

    // Save references to mesh and element base.
    mMesh = mesh;
    mReferenceQuad = quad;

    // Attach elements to mesh.
    mMesh->setupGlobalDof(mReferenceQuad->NumberDofVertex(),
                          mReferenceQuad->NumberDofEdge(),
                          mReferenceQuad->NumberDofFace(),
                          mReferenceQuad->NumberDofVolume(),
                          mReferenceQuad->NumberDimensions());

    // Setup boundary conditions from options.
    mMesh->setupBoundaries(options);

    // Register all global field required for time stepping.
    for (auto field : mMesh->GlobalFields()) {
        mMesh->registerFieldVectors(field);
    }

    // Get a list of all local elements.
    for (int i = 0; i < mMesh->NumberElementsLocal(); i++) {
        mElements.push_back(mReferenceQuad->clone());
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

        // Read auto-generated test function derivatives.
        element->readGradientOperator();

        // Assemble the (elemental) mass matrix.
        element->assembleElementMassMatrix(mMesh);

        // Attach any external source terms.
        element->attachSource(sources);

    }

    // Scatter the mass matrix to the global dofs.
    mMesh->checkInFieldBegin("m");
    mMesh->checkInFieldEnd("m");

}

void NewmarkGeneral::solve(Options options) {

    std::cout << "SOLVE!" << std::endl;

}
