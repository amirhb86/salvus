//
// Created by Michael Afanasiev on 2016-01-30.
//

#include <petscdm.h>
#include <iostream>
#include "SquareAcousticOrderFour.h"]
#include "Autogen/order4_square.h"

SquareAcousticOrderFour::SquareAcousticOrderFour(Options options) {

    // Basic properties.
    mElementShape = options.ElementShape();
    mPhysicsSystem = options.PhysicsSystem();
    mPolynomialOrder = options.PolynomialOrder();
    mNumberDimensions = 2;
    mNumberVertex = 4;

    // Gll points.
    mNumberDofVertex = 1;
    mNumberDofEdge = mPolynomialOrder - 1;
    mNumberDofFace = (mPolynomialOrder - 1) * (mPolynomialOrder - 1);
    mNumberDofVolume = 0;
    mNumberIntegrationPointsEps = 5;
    mNumberIntegrationPointsEta = 5;
    mNumberIntegrationPoints = 25;

    mIntegrationCoordinatesEps = {-1.0, -0.6546536707, 0.0, 0.6546536707, 1.0};
    mIntegrationCoordinatesEta = {-1.0, -0.6546536707, 0.0, 0.6546536707, 1.0};

}

void SquareAcousticOrderFour::registerFieldVectors() {

    PetscReal zero = 0;
    DMCreateLocalVector(mDistributedMesh, &mDisplacementLocal);
    DMCreateLocalVector(mDistributedMesh, &mAccelerationLocal);
    DMCreateLocalVector(mDistributedMesh, &mVelocityLocal);
    VecSet(mDisplacementLocal, zero);
    VecSet(mAccelerationLocal, zero);
    VecSet(mVelocityLocal, zero);

    DMCreateGlobalVector(mDistributedMesh, &mDisplacementGlobal);
    DMCreateGlobalVector(mDistributedMesh, &mAccelerationGlobal);
    DMCreateGlobalVector(mDistributedMesh, &mVelocityGlobal);
    VecSet(mDisplacementGlobal, zero);
    VecSet(mAccelerationGlobal, zero);
    VecSet(mVelocityGlobal, zero);

    PetscObjectSetName((PetscObject) mDisplacementGlobal, "displacement");
    PetscObjectSetName((PetscObject) mAccelerationGlobal, "acceleration");
    PetscObjectSetName((PetscObject) mVelocityGlobal, "velocity");

}

void SquareAcousticOrderFour::constructStiffnessMatrix() {

    double determinant_jacobian;
    Eigen::Matrix<double,2,2> jacobian, inverse_jacobian;
    for (auto eta: mIntegrationCoordinatesEta) {
        for (auto eps: mIntegrationCoordinatesEps) {

            // Get and invert Jacobian.
            jacobian = jacobianAtPoint(eps, eta);
            inverse_jacobian = jacobian.inverse();
            determinant_jacobian = jacobian.determinant();


        }
    }


}

void SquareAcousticOrderFour::interpolateMaterialProperties(ExodusModel &model) {

    // TODO. Test that this function results in a linear interpolation.

    mMaterialDensity.resize(mNumberIntegrationPoints);
    mMaterialVelocity.resize(mNumberIntegrationPoints);

    Eigen::Vector4d velocity_at_nodes;
    Eigen::Matrix<double,2,4> vertex_coordinates = VertexCoordinates();
    for (auto i = 0; i < mNumberVertex; i++) {
        velocity_at_nodes(i) = model.getMaterialParameterAtPoint({vertex_coordinates(0, i), vertex_coordinates(1, i)},
                                                                 "velocity");
    }

    int i = 0;
    for (auto eta: mIntegrationCoordinatesEta) {
        for (auto eps: mIntegrationCoordinatesEps) {
            Eigen::Vector4d coefficients = interpolateShapeFunctions(eps, eta);
            mMaterialVelocity(i) = coefficients.dot(velocity_at_nodes);
            i++;
        }

    }

}

void SquareAcousticOrderFour::readOperators() {

    double epsilon_0 = mIntegrationCoordinatesEps[0];
    double epsilon_1 = mIntegrationCoordinatesEps[1];
    double epsilon_2 = mIntegrationCoordinatesEps[2];
    double epsilon_3 = mIntegrationCoordinatesEps[3];
    double epsilon_4 = mIntegrationCoordinatesEps[4];

    double eta_0 = mIntegrationCoordinatesEta[0];
    double eta_1 = mIntegrationCoordinatesEta[1];
    double eta_2 = mIntegrationCoordinatesEta[2];
    double eta_3 = mIntegrationCoordinatesEta[3];
    double eta_4 = mIntegrationCoordinatesEta[4];

    int i = 0;
    double eta = eta_0;
    mGradientOperator.resize(mNumberIntegrationPointsEta, mNumberIntegrationPointsEps);
    Eigen::MatrixXd test(mNumberIntegrationPointsEta, mNumberIntegrationPointsEps);
    for (auto eps: mIntegrationCoordinatesEps) {
        interpolate_eps_derivative_order4_square(eps, epsilon_0, epsilon_1, epsilon_2, epsilon_3, epsilon_4,
                                                 eta, eta_0, eta_1, eta_2, eta_3, eta_4,
                                                 test.data());
        mGradientOperator.row(i) = test.col(0);
        i++;
    }
    std::cout << mGradientOperator << std::endl;
}
