#pragma once

#include <petscvec.h>
#include "Utilities/Options.h"
#include "Element/Simplex/Triangle.h"
#include "Model/ExodusModel.h"
#include "Utilities/Utilities.h"
#include "Mesh/Mesh.h"


class AcousticTri : public Triangle {

    Eigen::Vector3d mMaterialVelocityAtVertices;
    Eigen::MatrixXd mElementStrain;

    // precomputed element stiffness matrix (with velocities)
    Eigen::MatrixXd mElementStiffnessMatrix;
    
public:

    AcousticTri(Options options);

    AcousticTri *clone() const { return new AcousticTri(*this); }

    void computeSurfaceTerm();
    void assembleElementMassMatrix(Mesh *mesh);
    void interpolateMaterialProperties(ExodusModel *model);

    Eigen::MatrixXd computeSourceTerm(double time);
    Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement);

    void setInitialCondition(Mesh* mesh, Eigen::VectorXd& pts_x,Eigen::VectorXd& pts_z,
                             double L, double x0, double z0);
        
    Eigen::VectorXd exactSolution(Eigen::VectorXd& pts_x,Eigen::VectorXd& pts_z,
                                  double L, double x0, double z0, double time);

    std::vector<std::string> PullElementalFields() const { return {"u"}; }
    std::vector<std::string> PushElementalFields() const { return {"a"}; }

    Eigen::MatrixXd ElementStiffnessMatrix() { return mElementStiffnessMatrix; }

    void buildStiffnessMatrix();

    void prepareStiffness() { buildStiffnessMatrix(); }

    /**
     * Setup initial conditions for tests
     */
    void setupTest(Mesh* mesh, Options options);

    /**
     * Check exact solution against current displacement
     */
    double checkTest(Mesh* mesh, Options options, const Eigen::MatrixXd &displacement, double time);
    
};
