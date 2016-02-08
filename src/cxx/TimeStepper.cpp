//
// Created by Michael Afanasiev on 2016-02-08.
//

#include "TimeStepper.h"

void TimeStepper::registerFieldVector(std::string name, std::string component, Vec &field) {

    global_field_vector new_field;
    new_field.name = name;
    new_field.component = component;
    new_field.field = field;

}

void TimeStepper::registerMesh(DM &distributed_mesh, PetscSection &section) {

    mDistributedMesh = distributed_mesh;
    mSection = section;

}

void TimeStepper::takeTimeStep() {

}
