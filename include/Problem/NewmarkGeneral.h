#pragma once

// stl.
#include <vector>
#include <memory>

// parents.
#include <Receiver/Receiver.h>
#include <Problem/Problem.h>
#include <Mesh/Mesh.h>


class NewmarkGeneral: public Problem {

private:

    std::unique_ptr<Mesh> mMesh;
    std::shared_ptr<Element> mReferenceElem;
    std::vector<std::unique_ptr<Element>> mElements;
    std::vector<std::unique_ptr<Receiver>> mRecs;

public:

    ~NewmarkGeneral() {};
    virtual void solve(Options options);
    virtual void initialize(std::unique_ptr<Mesh> mesh, std::unique_ptr<ExodusModel> const &model,
                            std::unique_ptr<Options> const &options);

};
