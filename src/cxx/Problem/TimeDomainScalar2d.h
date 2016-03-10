//
// Created by Michael Afanasiev on 2016-02-23.
//

#ifndef SALVUS_TIMEDOMAINSCALAR2D_H
#define SALVUS_TIMEDOMAINSCALAR2D_H

#include "Problem.h"

class TimeDomainScalar2d : public Problem {

private:

    ScalarNewmark2D *mMesh;
    Acoustic *mReferenceQuad;
    std::vector<Acoustic *> mElements;

    double mSimulationDuration;
    double mTimeStep;


public:


    virtual void solve(Options options);
    virtual void initialize(Mesh *mesh, ExodusModel *model, Quad *quad, Options options);

};

#include "Problem.h"

#endif //SALVUS_TIMEDOMAINSCALAR2D_H
