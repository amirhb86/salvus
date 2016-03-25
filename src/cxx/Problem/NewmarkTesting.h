#pragma once

#include "Problem.h"

class NewmarkTesting: public Problem {

private:

    Mesh *mMesh;
    Element2D *mReferenceElem;
    std::vector<Element2D *> mElements;

public:

    void solve(Options options);
    void initialize(Mesh *mesh, ExodusModel *model, Element2D *elem, Options options);

    void applyDirichletBoundary(Mesh *mesh,std::string fieldname, double value);
    
};
