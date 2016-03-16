//
// Created by Michael Afanasiev on 2016-02-23.
//

#include "ScalarNewmark2D.h"

// Required fields for timestepping.

void ScalarNewmark2D::advanceField(double dt) {
    
    double pre_factor_acceleration = (1.0/2.0) * dt;
    double pre_factor_displacement = (1.0/2.0) * (dt * dt);

    VecAXPBYPCZ(mFields["v"].glb, pre_factor_acceleration, pre_factor_acceleration, 1.0,
                mFields["a"].glb, mFields["a_"].glb);

    VecAXPBYPCZ(mFields["u"].glb, dt, pre_factor_displacement, 1.0,
                mFields["v"].glb, mFields["a"].glb);

    VecCopy(mFields["a"].glb, mFields["a_"].glb);

}

void ScalarNewmark2D::applyInverseMassMatrix() {

    if (mFields.find("mi") == mFields.end()){
        registerFieldVectors("mi");
        VecCopy(mFields["m"].glb, mFields["mi"].glb);
        VecReciprocal(mFields["mi"].glb);
    }

    VecPointwiseMult(mFields["a"].glb, mFields["mi"].glb,
                     mFields["a"].glb);


}

std::vector<std::string> ScalarNewmark2D::GlobalFields() const {
    return {"u", "v", "a", "a_", "m"};
}