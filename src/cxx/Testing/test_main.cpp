//
// Created by Michael Afanasiev on 2016-03-27.
//

#define CATCH_CONFIG_RUNNER
#include "catch.h"
#include <Eigen/Dense>
#include <petsc.h>
#include <Element/Element.h>
#include <Element/HyperCube/Quad/Acoustic.h>

int main(int argc, char *argv[]) {

    // Init Salvus command line arguments.
    PetscInitialize(&argc, &argv, NULL, NULL);

    /* Need to be careful here. Catch defines its own command line
     * arguments, as does SALVUS. This is just a trick to get the
     * PETSc arguments read in, and then tell Catch that we really
     * don't have any command line arguments for it in particular.
     * */
    argc = 1;

    // Run all unit tests.
    int result = Catch::Session().run(argc, argv);

    // Clean up PETSc.
    PetscFinalize();

    return result;
}

TEST_CASE("Test whether simple stuff works.", "[element]") {

    // Common setup.
    Options options;
    options.setOptions();
    Mesh *mesh = Mesh::factory(options);
    mesh->read(options);
    ExodusModel *model = new ExodusModel(options);
    model->initializeParallel();
    std::vector<Source*> sources = Source::factory(options);
    Element2D *reference_element = Element2D::factory(options);

    // Define variables.
    Mesh *mMesh;
    Element2D *mRefElem;
    std::vector<Element2D*> mElements;

    // Problem.
    mMesh = mesh;
    mRefElem = reference_element;
    mMesh->setupGlobalDof(mRefElem->NumberDofVertex(),
                          mRefElem->NumberDofEdge(),
                          mRefElem->NumberDofFace(),
                          0,
                          mRefElem->NumberDimensions());
    mMesh->setupBoundaries(options);

    for (auto field : mMesh->GlobalFields()) {
        mMesh->registerFieldVectors(field);
    }

    for (int i = 0; i < mMesh->NumberElementsLocal(); i++) {
        mElements.push_back(mRefElem->clone());
    }

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

    }

    Eigen::MatrixXd tField = Eigen::MatrixXd::Zero(reference_element->NumberIntegrationPoints(), 1);
    Eigen::MatrixXd result = mElements[0]->computeStiffnessTerm(tField);


}

