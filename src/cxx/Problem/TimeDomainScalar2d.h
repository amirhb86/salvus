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

    /**
     * Sets values on the "dirichlet" boundary (specified in exodus
     * mesh file) to a specific value.     
     * @param mesh [in] Current mesh
     * @param fieldname [in] Set this fieldname to a specific value (e.g., "force")
     * @param value [in] Dirichlet value --- sets fieldname node on boundary to this value. Usually zero.
     */
    void applyDirichletBoundary(Mesh *mesh, std::string fieldname, double value);

    /**
     * Builds 1st-order absorbing boundary specified by "absorbing" in exodus file.
     * @param mesh [in] Current mesh
     */
    void absorbingBoundary(Mesh *mesh);
    
    virtual void solve(Options options);
    virtual void initialize(Mesh *mesh, ExodusModel *model, Quad *quad, Options options);

};

#include "Problem.h"

#endif //SALVUS_TIMEDOMAINSCALAR2D_H
