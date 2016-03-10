//
// Created by Michael Afanasiev on 2016-02-23.
//

#include "ElasticNewmark2D.h"

void ElasticNewmark2D::advanceField(double dt) {

    double pre_factor_acceleration = (1.0/2.0) * dt;
    double pre_factor_displacement = (1.0/2.0) * (dt * dt);

    VecAXPBYPCZ(mFields["velocity_x"].glb, pre_factor_acceleration, pre_factor_acceleration, 1.0,
                mFields["acceleration_x"].glb, mFields["acceleration_x_"].glb);
    VecAXPBYPCZ(mFields["velocity_z"].glb, pre_factor_acceleration, pre_factor_acceleration, 1.0,
                mFields["acceleration_z"].glb, mFields["acceleration_z_"].glb);

    VecAXPBYPCZ(mFields["displacement_x"].glb, dt, pre_factor_displacement, 1.0,
                mFields["velocity_x"].glb, mFields["acceleration_x"].glb);
    VecAXPBYPCZ(mFields["displacement_z"].glb, dt, pre_factor_displacement, 1.0,
                mFields["velocity_z"].glb, mFields["acceleration_z"].glb);

    VecCopy(mFields["acceleration_x"].glb, mFields["acceleration_x_"].glb);
    VecCopy(mFields["acceleration_z"].glb, mFields["acceleration_z_"].glb);

}

void ElasticNewmark2D::applyInverseMassMatrix() {

    if (mFields.find("mass_matrix_inverse") == mFields.end()) {
        registerFieldVectors("mass_matrix_inverse");
        VecCopy(mFields["mass_matrix"].glb, mFields["mass_matrix_inverse"].glb);
        VecReciprocal(mFields["mass_matrix_inverse"].glb);
    }

    VecPointwiseMult(mFields["acceleration_x"].glb, mFields["mass_matrix_inverse"].glb,
                     mFields["force_x"].glb);
    VecPointwiseMult(mFields["acceleration_z"].glb, mFields["mass_matrix_inverse"].glb,
                     mFields["force_z"].glb);

}
