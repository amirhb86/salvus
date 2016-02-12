#include <iostream>
#include <vector>
#include <mpi.h>
#include <petscsys.h>
#include "Mesh.h"

#include "Utilities.h"
#include "Element/HyperCube/Quad/Acoustic.h"
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
    Quad *reference_element = new Acoustic(options);

    mesh->setupGlobalDof(reference_element->NumberDofVertex(), reference_element->NumberDofEdge(),
                         reference_element->NumberDofFace(), reference_element->NumberDofVolume(),
                         reference_element->NumberDimensions());
    mesh->registerFields();



    // Clone a list of all local elements.
    std::vector<Quad *> elements;
    for (auto i = 0; i < mesh->NumberElementsLocal(); i++) { elements.push_back(reference_element->clone()); }

    // Now things that only local elements are allowed to do.
    int element_number = 0;
    for (auto &element: elements) {
        element->SetLocalElementNumber(element_number++);
        element->attachVertexCoordinates(mesh->DistributedMesh());
        element->attachSource(sources);
        element->interpolateMaterialProperties(model);
        element->readOperators();
        element->assembleMassMatrix();
        element->scatterMassMatrix(mesh);
    }

    mesh->checkInMassMatrix();
    mesh->setUpMovie(options.OutputMovieFile());

    double time = 0;
    double timestep = options.TimeStep();
    while (time < options.Duration()) {

        // Pull down fields from global dof.
        mesh->checkOutFields();
        mesh->zeroFields();

        // Compute element-wise terms.
        for (auto &element: elements) {
            element->SetTime(time);
            element->checkOutFields(mesh);
            element->computeStiffnessTerm();
            element->computeSourceTerm();
            element->computeSurfaceTerm();
            element->checkInField(mesh);
        }

        // Sum fields back into global dof.
        mesh->checkInFieldsBegin();
        mesh->checkInFieldsEnd();

        // Invert mass matrix and take time step.
        mesh->applyInverseMassMatrix();
        mesh->advanceField();
        mesh->saveFrame();

        time += timestep;
        if (!MPI::COMM_WORLD.Get_rank()) std::cout << time << std::endl;
//        break;
    }


    mesh->finalizeMovie();
    PetscFinalize();
}
