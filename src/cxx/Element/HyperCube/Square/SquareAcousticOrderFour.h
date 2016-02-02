//
// Created by Michael Afanasiev on 2016-01-30.
//

#ifndef SALVUS_SQUAREACOUSTICORDERFOUR_H
#define SALVUS_SQUAREACOUSTICORDERFOUR_H


#include <petscvec.h>
#include "../../Element.h"
#include "../../../Options.h"
#include "../Square.h"
#include "../../../Model/ExodusModel.h"

class SquareAcousticOrderFour: public Square {

    Vec mDisplacementLocal;
    Vec mAccelerationLocal;
    Vec mVelocityLocal;

    Vec mDisplacementGlobal;
    Vec mAccelerationGlobal;
    Vec mVelocityGlobal;

    Eigen::VectorXd mMaterialVelocity;
    Eigen::VectorXd mMaterialDensity;

    Eigen::MatrixXd mGradientOperator;

public:

    SquareAcousticOrderFour(Options options);

    virtual Element *clone() const { return new SquareAcousticOrderFour(*this); }

    virtual void readOperators();
    virtual void registerFieldVectors();
    virtual void constructStiffnessMatrix();
    virtual void interpolateMaterialProperties(ExodusModel &model);

};


#endif //SALVUS_SQUAREACOUSTICORDERFOUR_H
