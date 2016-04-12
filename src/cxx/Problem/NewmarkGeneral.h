//
// Created by Michael Afanasiev on 2016-03-14.
//

#ifndef SALVUS_NEWMARKGENERAL_H
#define SALVUS_NEWMARKGENERAL_H

#include "Problem.h"

class NewmarkGeneral: public Problem {

private:

    Mesh *mMesh;
    std::shared_ptr<Element> mReferenceElem;
    std::vector<std::shared_ptr<Element>> mElements;
    std::vector<std::shared_ptr<Receiver>> mRecs;

public:

    ~NewmarkGeneral() {};
    virtual void solve(Options options);
    virtual void initialize(Mesh *mesh, ExodusModel *model, std::shared_ptr<Element> elem, Options options);

};


#endif //SALVUS_NEWMARKGENERAL_H
