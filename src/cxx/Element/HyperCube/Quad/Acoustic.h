//
// Created by Michael Afanasiev on 2016-01-30.
//

#ifndef SALVUS_SQUAREACOUSTICORDERFOUR_H
#define SALVUS_SQUAREACOUSTICORDERFOUR_H


#include <petscvec.h>
#include "Utilities/Options.h"
#include "Element/HyperCube/Quad.h"
#include "Model/ExodusModel.h"
#include "Utilities/Utilities.h"
#include "Mesh/Mesh.h"


class AcousticQuad : public Quad {

    Eigen::Vector4d mMaterialVelocityAtVertices;
    Eigen::MatrixXd mElementStrain;

public:

    AcousticQuad(Options options);

    AcousticQuad *clone() const { return new AcousticQuad(*this); }

    void computeSurfaceTerm();
    void assembleElementMassMatrix(Mesh *mesh);
    void interpolateMaterialProperties(ExodusModel *model);

    Eigen::MatrixXd computeSourceTerm(double time);
    Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement);

    void setInitialCondition(Mesh* mesh, Eigen::VectorXd& pts_x, Eigen::VectorXd& pts_z);
    Eigen::VectorXd exactSolution(Eigen::VectorXd& pts_x,Eigen::VectorXd& pts_z,double time);

    std::vector<std::string> PullElementalFields() const { return {"u"}; }
    std::vector<std::string> PushElementalFields() const { return {"a"}; }

};


#endif //SALVUS_SQUAREACOUSTICORDERFOUR_H
