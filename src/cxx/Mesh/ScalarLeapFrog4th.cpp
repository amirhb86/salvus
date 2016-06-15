#include <Mesh/ScalarLeapFrog4th.h>
void ScalarLeapFrog4th::advanceField(double dt) {
    
    double pre_factor_acceleration = (1.0/2.0) * dt;
    double pre_factor_displacement = (1.0/2.0) * (dt * dt);

    // Goal
    // unp1[n] = 2*un[n] - unm1[n] - dt^2*an0[n] + dt^4/12.0*Minv[n]*b2[n]
    // Step1
    // unp1[n] = 2*un[n] - unm1[n]
    VecAXPBYPCZ(mFields["unp1"]->glb, // C*Z <-
                2.0, -1.0, 0.0, // A, B, C
                mFields["u"]->glb, // A*X
                mFields["unm1"]->glb // B*Y
                );
    // Step2
    // 4th order
    // unp1[n] += - dt^2*an0[n] + dt^4/12.0*an1[n]
    VecAXPBYPCZ(mFields["unp1"]->glb, // C*Z <-
                -dt*dt, dt*dt*dt*dt/12.0, 1.0, // A, B, C
                mFields["a0"]->glb, // A*X
                mFields["a1"]->glb // B*Y
                );
    
    // 2nd order
    // unp1[n] += - dt^2*an0[n] + dt^4/12.0*an1[n]
    // VecAXPBYPCZ(mFields["unp1"]->glb, // C*Z <-
    // -dt*dt, 0.0, 1.0, // A, B, C
    // mFields["a0"]->glb, // A*X
    // mFields["a1"]->glb // B*Y
    // );
    
    // unm1 <- un
    VecCopy(mFields["u"]->glb,mFields["unm1"]->glb);

    // un <- unp1
    VecCopy(mFields["unp1"]->glb,mFields["u"]->glb);
        
}

void ScalarLeapFrog4th::applyInverseMassMatrix() {

  if (mFields.find("mi") == mFields.end()){
      
    registerFieldVectors("mi");
    VecCopy(mFields["m"]->glb, mFields["mi"]->glb);
    VecReciprocal(mFields["mi"]->glb);
  }
    
  VecPointwiseMult(mFields["a"]->glb, mFields["mi"]->glb,
                   mFields["a"]->glb);

}

void ScalarLeapFrog4th::applyInverseMassMatrix(std::string tofield) {

  if (mFields.find("mi") == mFields.end()){
      
    registerFieldVectors("mi");
    VecCopy(mFields["m"]->glb, mFields["mi"]->glb);
    VecReciprocal(mFields["mi"]->glb);
  }
    
  VecPointwiseMult(mFields[tofield]->glb, mFields["mi"]->glb,
                   mFields[tofield]->glb);

}
