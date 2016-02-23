//
// Created by Michael Afanasiev on 2016-02-23.
//

#ifndef SALVUS_ELASTICNEWMARK2D_H
#define SALVUS_ELASTICNEWMARK2D_H

#include "Mesh.h"

class ElasticNewmark2D: public Mesh {

public:

    virtual void registerFields() {
        registerFieldVectors("acceleration_x_");
        registerFieldVectors("acceleration_z_");
        registerFieldVectors("acceleration_x");
        registerFieldVectors("acceleration_z");
        registerFieldVectors("displacement_x");
        registerFieldVectors("displacement_z");
        registerFieldVectors("velocity_x");
        registerFieldVectors("velocity_z");
        registerFieldVectors("force_x");
        registerFieldVectors("force_z");
        registerFieldVectors("mass_matrix");
    }

    virtual void advanceField();
    virtual void applyInverseMassMatrix();


};

#endif //SALVUS_ELASTICNEWMARK2D_H
