//
// Created by Michael Afanasiev on 2016-01-30.
//

#ifndef SALVUS_SQUAREACOUSTICORDERFOUR_H
#define SALVUS_SQUAREACOUSTICORDERFOUR_H


#include <petscvec.h>
#include "../../../Options.h"
#include "../Quad.h"
#include "../../../Model/ExodusModel.h"
#include "../../../Utilities.h"

class Acoustic : public Quad {

    Eigen::Vector4d mMaterialVelocityAtVertices;
    Eigen::Vector4d mMaterialDensityAtVertices;

    Eigen::VectorXd mIntegratedSource;
    Eigen::VectorXd mElementDisplacement;
    Eigen::VectorXd mIntegratedStiffnessMatrix;

    Eigen::MatrixXd mElementStress;
    Eigen::MatrixXd mElementStrain;

public:

    Acoustic(Options options);

    virtual Acoustic *clone() const { return new Acoustic(*this); }

    virtual void checkInField(Mesh *mesh);
    virtual void computeSourceTerm();
    virtual void computeSurfaceTerm();
    virtual void assembleMassMatrix();
    virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement);
    virtual void interpolateMaterialProperties(ExodusModel *model);

};


#endif //SALVUS_SQUAREACOUSTICORDERFOUR_H
