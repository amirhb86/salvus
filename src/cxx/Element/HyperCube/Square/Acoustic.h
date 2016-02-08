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

class Acoustic : public Square {

    Eigen::Vector4d mMaterialVelocityAtVertices;
    Eigen::Vector4d mMaterialDensityAtVertices;

    Eigen::VectorXd mElementForce;
    Eigen::VectorXd mElementDisplacement;
    Eigen::VectorXd mIntegratedStiffnessMatrix;
    Eigen::MatrixXd mElementStrain;

    double evaluateShapeFunctions(const double &eps, const double &eta, const int &itr);

public:

    Acoustic(Options options);

    virtual Element *clone() const { return new Acoustic(*this); }

    virtual void registerFieldVectors(Mesh *mesh);
    virtual void constructStiffnessMatrix(Mesh *mesh);
    virtual void interpolateMaterialProperties(ExodusModel &model);

};


#endif //SALVUS_SQUAREACOUSTICORDERFOUR_H
