//
// Created by Michael Afanasiev on 2016-02-08.
//

#ifndef SALVUS_TIMESTEPPER_H
#define SALVUS_TIMESTEPPER_H


#include <iosfwd>
#include <string>
#include <petscvec.h>
#include <vector>
#include <petscdmtypes.h>

class TimeStepper {

    struct global_field_vector {

        std::string name;
        std::string component;
        Vec field;
    };

    std::vector<global_field_vector> mFields;

    DM mDistributedMesh;
    PetscSection mSection;

    Vec mAccelerationScalar_;
    Vec mAccelerationX_;
    Vec mAccelerationY_;
    Vec mAccelerationZ_;

public:

    void registerMesh(DM &distributed_mesh, PetscSection &section);
    void registerFieldVector(std::string name, std::string component, Vec &field);
    void takeTimeStep();

};


#endif //SALVUS_TIMESTEPPER_H
