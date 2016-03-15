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
        registerFieldVectors("displacement_exact");
        registerFieldVectors("velocity");
        registerFieldVectors("force");
        registerFieldVectors("mass_matrix");
        registerFieldVectors("nodes_x");
        registerFieldVectors("nodes_z");
    }

    virtual void advanceField(double dt);
    virtual void applyInverseMassMatrix();
    virtual std::vector<std::string> GlobalFields() const;

};

#endif //SALVUS_SCALARNEWMARK2D_H
