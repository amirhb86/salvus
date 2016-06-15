#pragma once

#include <Mesh/Mesh.h>

class ScalarNewmark2D : public Mesh {

public:


    ~ScalarNewmark2D() {};
    
    ScalarNewmark2D(const std::unique_ptr<Options> &options): Mesh(options) { mGlobalFields = {"u", "v", "a", "a_", "m"}; }
//    ScalarNewmark2D(std::vector<std::string> fields) {mGlobalFields = fields;}
    
    virtual void advanceField(double dt);
    virtual void applyInverseMassMatrix();

};
