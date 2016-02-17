//
// Created by Michael Afanasiev on 2016-01-27.
//

#ifndef SALVUS_SOLVER_H
#define SALVUS_SOLVER_H


#include <ostream>
#include <iostream>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <iosfwd>
#include <string>
#include "Utilities.h"
#include "Element/HyperCube/Quad.h"
#include "Mesh.h"


class Problem {

public:

    static Problem *factory(std::string solver_type);

    virtual void solve() = 0;
    virtual void initialize(Mesh *mesh, ExodusModel *model, Quad *quad, Options options) = 0;


};

class TimeDomain : public Problem {

private:

    ScalarNewmark *mMesh;
    Quad *mReferenceQuad;
    std::vector<Quad *> mElements;

    double mSimulationDuration;
    double mTimeStep;


public:

    virtual void solve();
    virtual void initialize(Mesh *mesh, ExodusModel *model, Quad *quad, Options options);

};

class FrequencyDomain : public Problem {

public:

//    virtual void solve();
//    virtual void initialize(Mesh &mesh, ExodusModel &model, Quad &quad) = 0;

};


#endif //SALVUS_SOLVE
