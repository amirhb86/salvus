//
// Created by Michael Afanasiev on 2016-02-23.
//

#ifndef SALVUS_TIMEDOMAINELASTIC2D_H
#define SALVUS_TIMEDOMAINELASTIC2D_H

#include "Problem.h"

class TimeDomainElastic2d: public Problem {

private:

    ElasticNewmark2D *mMesh;
    Elastic *mReferenceQuad;
    std::vector<Elastic *> mElements;

    double mSimulationDuration;
    double mTimeStep;

public:

    virtual void solve();
    virtual void initialize(Mesh *mesh, ExodusModel *model, Quad *quad, Options options);


};

#include "Problem.h"

#endif //SALVUS_TIMEDOMAINELASTIC2D_H
