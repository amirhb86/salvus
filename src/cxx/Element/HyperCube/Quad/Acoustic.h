//
// Created by Michael Afanasiev on 2016-01-30.
//

#ifndef SALVUS_SQUAREACOUSTICORDERFOUR_H
#define SALVUS_SQUAREACOUSTICORDERFOUR_H


#include <petscvec.h>
#include "../../../Utilities/Options.h"
#include "../Quad.h"
#include "../../../Model/ExodusModel.h"
#include "../../../Utilities/Utilities.h"
#include "../../../Mesh/Mesh.h"

class Acoustic : public Quad {

    Eigen::Vector4d mMaterialVelocityAtVertices;
    Eigen::Vector4d mMaterialDensityAtVertices;
    Eigen::MatrixXd mElementStrain;

public:

    Acoustic(Options options);

    virtual Acoustic *clone() const { return new Acoustic(*this); }

    virtual void computeSurfaceTerm();
    virtual void assembleElementMassMatrix(Mesh *mesh);
    virtual void interpolateMaterialProperties(ExodusModel *model);

    virtual Eigen::MatrixXd computeSourceTerm();
    virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement);

    void setInitialCondition(Mesh* mesh, Eigen::VectorXd& pts_x, Eigen::VectorXd& pts_z);
    Eigen::VectorXd exactSolution(Eigen::VectorXd& pts_x,Eigen::VectorXd& pts_z,double time);
};


#endif //SALVUS_SQUAREACOUSTICORDERFOUR_H
