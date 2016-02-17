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
    virtual void checkOutField(Mesh *mesh);
    virtual Eigen::MatrixXd computeSourceTerm();
    virtual void computeSurfaceTerm();
    virtual void assembleElementMassMatrix(Mesh *mesh);
    virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement);
    virtual void interpolateMaterialProperties(ExodusModel *model);

    virtual void setInitialCondition(Mesh* mesh, Eigen::VectorXd& pts_x, Eigen::VectorXd& pts_z);
    
};


#endif //SALVUS_SQUAREACOUSTICORDERFOUR_H
