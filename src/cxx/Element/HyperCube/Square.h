//
// Created by Michael Afanasiev on 2016-01-30.
//

#ifndef SALVUS_HYPERCUBE_H
#define SALVUS_HYPERCUBE_H

#include "../Element.h"
#include "../../Model/ExodusModel.h"

extern "C" {
#include "Square/Autogen/order4_square.h"
};

/*
 * Base class of an abstract four node square. The reference element is set up as below.
 *
 * (n2)______________(n3)
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * (n0)______________(n1)
 *
 * (eta)
 *   ^
 *   |
 *   |
 *   |______> (eps)
*/


class Square : public Element {

    static PetscReal n0(const PetscReal eps, const PetscReal eta);
    static PetscReal n1(const PetscReal eps, const PetscReal eta);
    static PetscReal n2(const PetscReal eps, const PetscReal eta);
    static PetscReal n3(const PetscReal eps, const PetscReal eta);

    static PetscReal dn0deps(const PetscReal eta);
    static PetscReal dn1deps(const PetscReal eta);
    static PetscReal dn2deps(const PetscReal eta);
    static PetscReal dn3deps(const PetscReal eta);

    static PetscReal dn0deta(const PetscReal eps);
    static PetscReal dn1deta(const PetscReal eps);
    static PetscReal dn2deta(const PetscReal eps);
    static PetscReal dn3deta(const PetscReal eps);

    std::vector<PetscReal> mIntegrationPoints;

    // Fixed size eigen matrices for speed
    bool mCheckHull(double x, double z);
    Eigen::Matrix<double,2,4> mVertexCoordinates;

    Eigen::Vector2d inverseCoordinateTransform(const double &x_real, const double &z_real,
                                               double eps, double eta);

protected:

    PetscInt mNumberIntegrationPointsEps;
    PetscInt mNumberIntegrationPointsEta;

    Eigen::VectorXd mIntegrationCoordinatesEps;
    Eigen::VectorXd mIntegrationCoordinatesEta;
    Eigen::VectorXd mIntegrationWeightsEps;
    Eigen::VectorXd mIntegrationWeightsEta;

    Eigen::MatrixXd mGradientOperator;

public:

    // Local methods.
    Eigen::Vector4d __interpolateMaterialProperties(ExodusModel &model, std::string parameter_name);
    Eigen::Matrix<double,2,2> jacobianAtPoint(PetscReal eps, PetscReal eta);
    Eigen::Vector4d interpolateShapeFunctions(PetscReal eps, PetscReal eta);

    Eigen::Map<Eigen::VectorXd> epsVectorStride(
            Eigen::VectorXd &function, int &eta_index);
    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<>> etaVectorStride(
            Eigen::VectorXd &function, int &eta_index);

    // Utility methods.
    virtual Element* clone() const = 0;

    // Superclass methods
    virtual void readOperators();
    virtual void attachVertexCoordinates();
    virtual void attachIntegrationPoints();
    virtual void attachSource(std::vector<Source*> sources);

    // Pure virtual methods.
    virtual void registerFieldVectors(Mesh *mesh) = 0;
    virtual void constructStiffnessMatrix(Mesh *mesh) = 0;
    virtual void interpolateMaterialProperties(ExodusModel &model) = 0;

    // Getattrs.
    Eigen::Matrix<double,2,4> VertexCoordinates() { return mVertexCoordinates; }

};


#endif //SALVUS_HYPERCUBE_H
