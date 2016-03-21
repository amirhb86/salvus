//
// Created by Michael Afanasiev on 2016-02-23.
//

#ifndef SALVUS_ELASTICNEWMARK2D_H
#define SALVUS_ELASTICNEWMARK2D_H

#include "Mesh.h"

class ElasticNewmark2D: public Mesh {

public:

    virtual void advanceField(double dt);
    virtual void applyInverseMassMatrix();
    virtual std::vector<std::string> GlobalFields() const;


};

#endif //SALVUS_ELASTICNEWMARK2D_H
