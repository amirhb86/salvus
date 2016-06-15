#pragma once

#include <Mesh/Mesh.h>

class ElasticNewmark2D: public Mesh {

public:

    ~ElasticNewmark2D() {};
    virtual void advanceField(double dt);
    virtual void applyInverseMassMatrix();

    ElasticNewmark2D(const std::unique_ptr<Options> &options): Mesh(options) {
      mGlobalFields = {"ux", "uy", "vx", "vy", "ax", "ay", "ax_", "ay_", "m"};
    }
//    ElasticNewmark2D(std::vector<std::string> fields) {mGlobalFields = fields;}
    
};
