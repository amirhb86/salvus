//
// Created by Michael Afanasiev on 2016-02-23.
//

#include "TimeDomainElastic2d.h"

void TimeDomainElastic2d::initialize(Mesh *mesh, ExodusModel *model, Quad *quad, Options options) {

    // Get sources.
    std::vector<Source*> sources = Source::factory(options);

    // Perform dynamic casts to ensure types are appropriate.
    mMesh = dynamic_cast<ElasticNewmark2D *> (mesh);
    mReferenceQuad = dynamic_cast<Elastic*> (quad);


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
        element->attachSource(sources);
        element->readGradientOperator();
        element->assembleElementMassMatrix(mesh);
    }

    // Assemble the mass matrix to the global dof.
    mesh->assembleLocalFieldToGlobal("mass_matrix");

    // Set up options.
    mMesh->setUpMovie(options.OutputMovieFile());
    mSimulationDuration = options.Duration();
    mTimeStep = options.TimeStep();

}

void TimeDomainElastic2d::solve(Options options) {


    Eigen::MatrixXd F;
    Eigen::MatrixXd Ku;
    Eigen::MatrixXd FminusKu;
    Eigen::MatrixXd displacement_on_element(2, mReferenceQuad->NumberIntegrationPoints()); // x-displacement -> row(0), z-displacement -> row(1)
    double time = 0;
    int it = 0;
    while (time < mSimulationDuration) {

        mMesh->checkOutField("displacement_x");
        mMesh->checkOutField("displacement_z");
        mMesh->zeroFields("force_x");
        mMesh->zeroFields("force_z");

        for (auto &element: mElements) {

            displacement_on_element.row(0) = element->checkOutFieldElement(mMesh, "displacement_x");
            displacement_on_element.row(1) = element->checkOutFieldElement(mMesh, "displacement_z");

            F = element->computeSourceTerm(time);
            Ku = element->computeStiffnessTerm(displacement_on_element);
            FminusKu = F - Ku;

            element->checkInFieldElement(mMesh, FminusKu.row(0), "force_x");
            element->checkInFieldElement(mMesh, FminusKu.row(1), "force_z");

        }

        mMesh->checkInFieldBegin("force_x");
        mMesh->checkInFieldBegin("force_z");
        mMesh->checkInFieldEnd("force_x");
        mMesh->checkInFieldEnd("force_z");

        mMesh->applyInverseMassMatrix();
        mMesh->advanceField(mTimeStep);

        mMesh->saveFrame("displacement_x",it);

        std::cout << time << std::endl;
        time += mTimeStep;
        it++;
    }

}
