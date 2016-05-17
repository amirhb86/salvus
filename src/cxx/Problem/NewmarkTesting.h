#pragma once

#include "Problem.h"

class NewmarkTesting: public Problem {

private:

    Mesh *mMesh;
    std::shared_ptr<ElementNew> mReferenceElem;
    std::vector<std::shared_ptr<ElementNew>> mElements;

public:

    void solve(Options options);
    void initialize(Mesh *mesh,
                    ExodusModel *model,
                    std::shared_ptr<ElementNew> elem, Options options);

    
};
