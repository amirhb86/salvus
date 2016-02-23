//
// Created by Michael Afanasiev on 2016-02-23.
//

#include "ScalarNewmark2D.h"

void ScalarNewmark2D::advanceField() {

    double dt = 1e-3;
    double pre_factor_acceleration = (1.0/2.0) * dt;
    double pre_factor_displacement = (1.0/2.0) * (dt * dt);

    VecAXPBYPCZ(mFields["velocity"].glb, pre_factor_acceleration, pre_factor_acceleration, 1.0,
                mFields["acceleration"].glb, mFields["acceleration_"].glb);

    VecAXPBYPCZ(mFields["displacement"].glb, dt, pre_factor_displacement, 1.0,
                mFields["velocity"].glb, mFields["acceleration"].glb);

    VecCopy(mFields["acceleration"].glb, mFields["acceleration_"].glb);

}

void ScalarNewmark2D::applyInverseMassMatrix() {

    if (mFields.find("mass_matrix_inverse") == mFields.end()){
        registerFieldVectors("mass_matrix_inverse");
        VecCopy(mFields["mass_matrix"].glb, mFields["mass_matrix_inverse"].glb);
        VecReciprocal(mFields["mass_matrix_inverse"].glb);
    }

    VecPointwiseMult(mFields["acceleration"].glb, mFields["mass_matrix_inverse"].glb,
                     mFields["force"].glb);

}