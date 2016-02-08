#include <iostream>
#include <vector>
#include <mpi.h>
#include <petscsys.h>
#include "Mesh.h"

#include "Element/Element.h"
#include "Utilities.h"
#include "Element/HyperCube/Square/Acoustic.h"
#include "Model/ExodusModel.h"
#include "Source.h"

static constexpr char help[] = "Welcome to salvus.";

int main(int argc, char *argv[]) {

    PetscInitialize(&argc, &argv, NULL, help);

    // Get command line options.
    Options options;
    options.setOptions();

    // Get mesh.
    Mesh *mesh = Mesh::factory(options);
    mesh->read(options);

    // Get model.
    ExodusModel model(options);
    model.initializeParallel();

    // Get timestepper.
    TimeStepper *time_stepper = new TimeStepper;

    // Get source.
    std::vector<Source*> sources = Source::factory(options);

    // Setup reference element.
    Element *reference_element = Element::factory(options);
    mesh->setupGlobalDof(reference_element->NumberDofVertex(),
                         reference_element->NumberDofEdge(),
                         reference_element->NumberDofFace(),
                         reference_element->NumberDofVolume(),
                         reference_element->NumberDimensions());

    time_stepper->registerMesh(mesh->DistributedMesh(), mesh->MeshSection());

    // Now that the mesh is constructed, register it with the reference element.
    reference_element->registerMesh(mesh->DistributedMesh(), mesh->MeshSection());
    reference_element->registerFieldVectors(mesh);

    // Clone a list of all local elements.
    std::vector<Element *> elements;
    for (auto i = 0; i < mesh->NumberElementsLocal(); i++) {
        elements.push_back(reference_element->clone());
    }

    // Now things that only local elements are allowed to do.
    PetscInt element_number = 0;
    for (auto &element: elements) {
        element->SetLocalElementNumber(element_number++);
        element->attachVertexCoordinates();
        element->attachIntegrationPoints();
        element->attachSource(sources);
        element->interpolateMaterialProperties(model);
        element->readOperators();
    }


    while (true) {
        mesh->checkOutFields();
        for (auto &element: elements) { element->constructStiffnessMatrix(mesh); }
        mesh->checkInFieldsBegin();
        mesh->checkInFieldsEnd();
        mesh->advanceField();
        break;
    }


    PetscFinalize();
}
