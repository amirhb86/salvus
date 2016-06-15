#include <Mesh/ElasticNewmark2D.h>

void ElasticNewmark2D::advanceField(double dt) {

    double pre_factor_acceleration = (1.0/2.0) * dt;
    double pre_factor_displacement = (1.0/2.0) * (dt * dt);

    VecAXPBYPCZ(mFields["vx"]->glb, pre_factor_acceleration, pre_factor_acceleration, 1.0,
                mFields["ax"]->glb, mFields["ax_"]->glb);
    VecAXPBYPCZ(mFields["vy"]->glb, pre_factor_acceleration, pre_factor_acceleration, 1.0,
                mFields["ay"]->glb, mFields["ay_"]->glb);

    VecAXPBYPCZ(mFields["ux"]->glb, dt, pre_factor_displacement, 1.0,
                mFields["vx"]->glb, mFields["ax"]->glb);
    VecAXPBYPCZ(mFields["uy"]->glb, dt, pre_factor_displacement, 1.0,
                mFields["vy"]->glb, mFields["ay"]->glb);

    VecCopy(mFields["ax"]->glb, mFields["ax_"]->glb);
    VecCopy(mFields["ay"]->glb, mFields["ay_"]->glb);

}

void ElasticNewmark2D::applyInverseMassMatrix() {

    if (mFields.find("mi") == mFields.end()) {
        registerFieldVectors("mi");
        VecCopy(mFields["m"]->glb, mFields["mi"]->glb);
        VecReciprocal(mFields["mi"]->glb);
    }

    VecPointwiseMult(mFields["ax"]->glb, mFields["mi"]->glb,
                     mFields["ax"]->glb);
    VecPointwiseMult(mFields["ay"]->glb, mFields["mi"]->glb,
                     mFields["ay"]->glb);

}

