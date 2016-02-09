#include <iostream>
#include <vector>
#include <mpi.h>
#include <petscsys.h>
#include "Mesh.h"

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

    // Get source.
    std::vector<Source*> sources = Source::factory(options);

    // Setup reference element.
    Square *reference_element = new Acoustic(options);
    mesh->setupGlobalDof(reference_element->NumberDofVertex(), reference_element->NumberDofEdge(),
                         reference_element->NumberDofFace(), reference_element->NumberDofVolume(),
                         reference_element->NumberDimensions());


    // Clone a list of all local elements.
    std::vector<Square *> elements;
    for (auto i = 0; i < mesh->NumberElementsLocal(); i++) { elements.push_back(reference_element->clone()); }

    // Now things that only local elements are allowed to do.
    int element_number = 0;
    for (auto &element: elements) {
        element->SetLocalElementNumber(element_number++);
        element->attachVertexCoordinates(mesh->DistributedMesh());
        element->attachSource(sources);
        element->interpolateMaterialProperties(model);
        element->readOperators();
    }


    mesh->registerFields();
    while (true) {
        mesh->checkOutFields();
        for (auto &element: elements) {
            element->checkOutFields(mesh);
            element->computeStiffnessTerm();
            element->computeSourceTerm();
            element->computeSurfaceTerm();
            element->checkInField(mesh);
        }
        mesh->checkInFieldsBegin();
        mesh->checkInFieldsEnd();
        mesh->advanceField();
        std::cout << "LOOP" << std::endl;
        break;
    }


    PetscFinalize();
}
