//
// Created by Michael Afanasiev on 2016-02-23.
//

#ifndef SALVUS_SCALARNEWMARK2D_H
#define SALVUS_SCALARNEWMARK2D_H

#include "Mesh.h"

class ScalarNewmark2D : public Mesh {

public:

    virtual void registerFields() {
        registerFieldVectors("acceleration_");
        registerFieldVectors("acceleration");
        registerFieldVectors("displacement");
        registerFieldVectors("velocity");
        registerFieldVectors("force");
        registerFieldVectors("mass_matrix");
        registerFieldVectors("nodes_x");
        registerFieldVectors("nodes_z");
    }

    virtual void advanceField();
    virtual void applyInverseMassMatrix();

};

#endif //SALVUS_SCALARNEWMARK2D_H
