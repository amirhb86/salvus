#pragma once

#include "Problem.h"

class NewmarkTesting: public Problem {

private:

    Mesh *mMesh;
    Element *mReferenceElem;
    std::vector<Element *> mElements;

public:

    void solve(Options options);
    void initialize(Mesh *mesh, ExodusModel *model, Element *elem, Options options);

    
};
