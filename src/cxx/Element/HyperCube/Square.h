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
    Eigen::Matrix<double,2,4> mVertexCoordinates;

protected:

    PetscInt mNumberIntegrationPointsEps;
    PetscInt mNumberIntegrationPointsEta;

    std::vector<PetscReal> mIntegrationCoordinatesEps;
    std::vector<PetscReal> mIntegrationCoordinatesEta;

    Eigen::VectorXd mIntegrationWeightsEps;
    Eigen::VectorXd mIntegrationWeightsEta;
    Eigen::MatrixXd mGradientOperator;

public:

    // Local methods.
    Eigen::Matrix<double,2,2> jacobianAtPoint(PetscReal eps, PetscReal eta);
    Eigen::Vector4d interpolateShapeFunctions(PetscReal eps, PetscReal eta);

    Eigen::Map<Eigen::VectorXd> epsVectorStride(
            Eigen::VectorXd &function, int &eta_index);
    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<>> etaVectorStride(
            Eigen::VectorXd &function, int &eta_index);

    // Utility methods.
    virtual Element* clone() const = 0;

    // Superclass methods
    virtual void attachVertexCoordinates();
    virtual void attachIntegrationPoints();

    // Pure virtual methods.
    virtual void readOperators() = 0;
    virtual void registerFieldVectors() = 0;
    virtual void constructStiffnessMatrix() = 0;
    virtual void interpolateMaterialProperties(ExodusModel &model) = 0;

    // Getattrs.
    Eigen::Matrix<double,2,4> VertexCoordinates() { return mVertexCoordinates; }

};


#endif //SALVUS_HYPERCUBE_H
