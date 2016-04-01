//
// Created by Michael Afanasiev on 2016-01-30.
//

#pragma once

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
    void attachMaterialProperties(ExodusModel *model);

    Eigen::MatrixXd computeSourceTerm(double time);
    Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement);

    void setInitialCondition(Mesh* mesh, Eigen::VectorXd& pts_x,Eigen::VectorXd& pts_z,
                             double L, double x0, double z0);
        
    Eigen::VectorXd exactSolution(Eigen::VectorXd& pts_x,Eigen::VectorXd& pts_z,
                                  double L, double x0, double z0, double time);

    std::vector<std::string> PullElementalFields() const { return {"u"}; }
    std::vector<std::string> PushElementalFields() const { return {"a"}; }

    /**
     * Setup initial conditions for tests
     */
    void setupTest(Mesh* mesh, Options options);

    /**
     * Check exact solution against current displacement
     */
    double checkTest(Mesh* mesh, Options options, const Eigen::MatrixXd &displacement, double time);
    
    /**
     * Empty as its not needed for quadrilaterals.
     */ 
    void prepareStiffness() {  }
    
};

