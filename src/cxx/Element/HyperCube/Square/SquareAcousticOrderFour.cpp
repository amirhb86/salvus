//
// Created by Michael Afanasiev on 2016-01-30.
//

#include <petscdm.h>
#include <iostream>
#include <petscdmplex.h>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include "SquareAcousticOrderFour.h"
#include "Autogen/order4_square.h"

SquareAcousticOrderFour::SquareAcousticOrderFour(Options options) {


    // Basic properties.
    SetNumberVertex(4);
    SetNumberDimensions(2);
    SetElementShape(options.ElementShape());
    SetPhysicsSystem(options.PhysicsSystem());
    SetPolynomialOrder(options.PolynomialOrder());
    SetContainsSource(false);

    // Gll points.
    mNumberDofVertex = 1;
    mNumberDofEdge = mPolynomialOrder - 1;
    mNumberDofFace = (mPolynomialOrder - 1) * (mPolynomialOrder - 1);
    mNumberDofVolume = 0;
    mNumberIntegrationPointsEps = 5;
    mNumberIntegrationPointsEta = 5;
    mNumberIntegrationPoints = 25;

    // Integration points.
    mIntegrationCoordinatesEps = {-1.0, -0.6546536707, 0.0, 0.6546536707, 1.0};
    mIntegrationCoordinatesEta = {-1.0, -0.6546536707, 0.0, 0.6546536707, 1.0};

    mIntegrationWeightsEta.resize(mNumberIntegrationPointsEta);
    mIntegrationWeightsEps.resize(mNumberIntegrationPointsEps);
    mIntegrationWeightsEps << 0.1, 0.5444444444, 0.7111111111, 0.5444444444, 0.1;
    mIntegrationWeightsEta << 0.1, 0.5444444444, 0.7111111111, 0.5444444444, 0.1;

    // Allocate element matrices.
    mElementForce.setZero(mNumberIntegrationPoints);
    mElementDisplacement.setZero(mNumberIntegrationPoints);
    mIntegratedStiffnessMatrix.setZero(mNumberIntegrationPoints);

    // Parameter mapping.
    mClosureMapping = {
            6, 13, 22, 3, 15,
            7, 16, 23, 2, 20,
            8, 17, 19, 1, 24,
            11, 18, 14, 5, 4,
            12, 21, 9, 10, 0};

}

void SquareAcousticOrderFour::registerFieldVectors() {

    PetscReal zero = 0;
    DMCreateLocalVector(mDistributedMesh, &mElementalForceLocal);
    DMCreateLocalVector(mDistributedMesh, &mDisplacementLocal);
    DMCreateLocalVector(mDistributedMesh, &mAccelerationLocal);
    DMCreateLocalVector(mDistributedMesh, &mVelocityLocal);
    VecSet(mDisplacementLocal, zero);
    VecSet(mAccelerationLocal, zero);
    VecSet(mVelocityLocal, zero);

    DMCreateGlobalVector(mDistributedMesh, &mElementalForceGlobal);
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

    // Test
    for (int i = 0; i < mNumberIntegrationPoints; i++) { mElementDisplacement(i) = mIntegrationCoordinatesEta[i%5]; }

    int itr = 0;
    Eigen::VectorXd divergence(mNumberIntegrationPoints);
    Eigen::VectorXd external_forcing(mNumberIntegrationPoints);
    Eigen::Matrix<double,2,1> strain;
    Eigen::Matrix<double,2,1> test_function_gradient;
    Eigen::Matrix<double,2,2> inverse_Jacobian;

    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            // Eps and eta coordinates.
            double eta = mIntegrationCoordinatesEta[eta_index];
            double eps = mIntegrationCoordinatesEps[eps_index];

            // Get and invert Jacobian.
            inverse_Jacobian = jacobianAtPoint(eps, eta).inverse();

            // Calculate strain.
            strain(0) = mGradientOperator.row(eps_index).dot(epsVectorStride(mElementDisplacement, eta_index));
            strain(1) = mGradientOperator.row(eta_index).dot(etaVectorStride(mElementDisplacement, eta_index));
            strain = inverse_Jacobian * strain;

            // Calculate test function derivatives.
            test_function_gradient(0) = mGradientOperator.row(eps_index).sum();
            test_function_gradient(1) = mGradientOperator.row(eta_index).sum();
            test_function_gradient = inverse_Jacobian * test_function_gradient;

            // Get material parameters at this node.
            double velocity = interpolateShapeFunctions(eps, eta).dot(mMaterialVelocityAtVertices);
            divergence(itr) = velocity * test_function_gradient.dot(strain);
            itr++;

        }
    }

    itr = 0;
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            double eps = mIntegrationCoordinatesEps[eps_index];
            double eta = mIntegrationCoordinatesEps[eta_index];
            double determinant = jacobianAtPoint(eps, eta).determinant();

            mIntegratedStiffnessMatrix(itr) = determinant *
                    (Square::mIntegrationWeightsEps.dot(epsVectorStride(divergence, eps_index)) +
                     mIntegrationWeightsEta.dot(etaVectorStride(divergence, eps_index)));

            for (auto &source: mSources) {
                external_forcing(itr) += evaluateShapeFunctions(
                        source->ReferenceLocationEps(), source->ReferenceLocationEta(), itr) /
                                         (Square::mIntegrationWeightsEps(eps_index) * mIntegrationWeightsEta(eta_index) *
                        determinant) * source->fire(mTime);
            }

            itr++;
        }
    }

    mElementForce = external_forcing - mIntegratedStiffnessMatrix;

}

void SquareAcousticOrderFour::interpolateMaterialProperties(ExodusModel &model) {

    mMaterialVelocityAtVertices = __interpolateMaterialProperties(model, "velocity");

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
}

void SquareAcousticOrderFour::gatherPartitionFieldsToElement() {
    __gatherPartitionFieldsToElement(mDisplacementLocal, mElementDisplacement, mClosureMapping);
}

void SquareAcousticOrderFour::scatterElementFieldsToPartition() {
    __scatterElementFieldsToPartition(mElementalForceLocal, mElementForce, mClosureMapping);
}

void SquareAcousticOrderFour::gatherDistributedFieldsToPartition() {
    __gatherDistributedFieldsToPartition(mDisplacementLocal, mDisplacementGlobal);
}

void SquareAcousticOrderFour::scatterPartitionFieldsToDistributedBegin() {
    __scatterPartitionFieldsToDistributedBegin(mElementalForceLocal, mElementalForceGlobal);
}

void SquareAcousticOrderFour::scatterPartitionFieldsToDistributedEnd() {
    __scatterPartitionFieldsToDistributedEnd(mElementalForceLocal, mElementalForceGlobal);
}

double SquareAcousticOrderFour::evaluateShapeFunctions(const double &eps, const double &eta, const int &itr) {

    Eigen::VectorXd test(mNumberIntegrationPoints);
    interpolate_order4_square(eps, mIntegrationCoordinatesEps[0],
                              mIntegrationCoordinatesEps[1], mIntegrationCoordinatesEps[2],
                              mIntegrationCoordinatesEps[3], mIntegrationCoordinatesEps[4], eta,
                              mIntegrationCoordinatesEta[0], mIntegrationCoordinatesEta[1],
                              mIntegrationCoordinatesEta[2], mIntegrationCoordinatesEta[3],
                              mIntegrationCoordinatesEta[4], test.data());

    return test(itr);

}

