#include <iostream>
#include <vector>
#include <mpi.h>
#include <petscsys.h>
#include "Mesh.h"

#include "Utilities.h"
#include "Element/HyperCube/Quad/Acoustic.h"
#include "Model/ExodusModel.h"
#include "Source.h"
#include "Problem.h"

static constexpr char help[] = "Welcome to salvus.";

int main(int argc, char *argv[]) {

    PetscInitialize(&argc, &argv, NULL, help);

    // Get command line options. DUMB THING IN NEED TO DO. OTHER DUMB THINGS.
    Options options;
    options.setOptions();

    // Get mesh.
    Mesh *mesh = Mesh::factory(options);
    mesh->read(options);

    // Get model.
    ExodusModel *model = new ExodusModel(options);
    model->initializeParallel();

    // Get sources.
    std::vector<Source*> sources = Source::factory(options);

    // Setup reference element.
    Quad *reference_element = new Acoustic(options);

    // Use above elements to define the problem.
    Problem *problem = Problem::factory("time_domain");
    problem->initialize(mesh, model, reference_element);
    problem->solve();

    PetscFinalize();
}
