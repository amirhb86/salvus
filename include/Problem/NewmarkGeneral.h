#pragma once

// stl.
#include <vector>
#include <memory>

// parents.
#include <Problem/Problem.h>

class Receiver;

class NewmarkGeneral: public Problem {

private:

    Mesh *mMesh;
    std::shared_ptr<Element> mReferenceElem;
    std::vector<std::shared_ptr<Element>> mElements;
    std::vector<std::shared_ptr<Receiver>> mRecs;

public:

    ~NewmarkGeneral() {};
    virtual void solve(Options options);
    virtual void initialize(Mesh *mesh, ExodusModel *model, std::unique_ptr<Options> const &options);

};
