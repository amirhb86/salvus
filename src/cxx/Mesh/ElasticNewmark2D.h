//
// Created by Michael Afanasiev on 2016-02-23.
//

#ifndef SALVUS_ELASTICNEWMARK2D_H
#define SALVUS_ELASTICNEWMARK2D_H

#include "Mesh.h"

class ElasticNewmark2D: public Mesh {

public:

    ~ElasticNewmark2D() {};
    virtual void advanceField(double dt);
    virtual void applyInverseMassMatrix();

    ElasticNewmark2D() { mGlobalFields = {"ux", "uy", "vx", "vy", "ax", "ay", "ax_", "ay_", "m"}; }
    ElasticNewmark2D(std::vector<std::string> fields) {mGlobalFields = fields;}
    
};

#endif //SALVUS_ELASTICNEWMARK2D_H
