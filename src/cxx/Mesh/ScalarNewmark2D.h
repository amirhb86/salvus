//
// Created by Michael Afanasiev on 2016-02-23.
//

#ifndef SALVUS_SCALARNEWMARK2D_H
#define SALVUS_SCALARNEWMARK2D_H

#include "Mesh.h"

class ScalarNewmark2D : public Mesh {

public:
    
    ScalarNewmark2D() { mGlobalFields = {"u", "v", "a", "a_", "m"}; }
    ScalarNewmark2D(std::vector<std::string> fields) {mGlobalFields = fields;}
    
    virtual void advanceField(double dt);
    virtual void applyInverseMassMatrix();

};

#endif //SALVUS_SCALARNEWMARK2D_H
