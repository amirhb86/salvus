//
// Created by Michael Afanasiev on 2016-01-30.
//

#include <petscdm.h>
#include <iostream>
#include <petscdmplex.h>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include "Acoustic.h"
#include "Autogen/order4_square.h"

Acoustic::Acoustic(Options options) {


    // Basic properties.
    SetNumberVertex(4);
    SetNumberDimensions(2);
    SetElementShape(options.ElementShape());
    SetPhysicsSystem(options.PhysicsSystem());
    SetPolynomialOrder(options.PolynomialOrder());
    SetContainsSource(false);

    // Gll points.
    mNumberDofVolume = 0;
    mNumberDofVertex = 1;
    mNumberDofEdge = mPolynomialOrder - 1;
    mNumberDofFace = (mPolynomialOrder - 1) * (mPolynomialOrder - 1);

    // Integration points.
    mIntegrationCoordinatesEps = Element::GllPointsForOrder(options.PolynomialOrder());
    mIntegrationCoordinatesEta = Element::GllPointsForOrder(options.PolynomialOrder());
    mIntegrationWeightsEps = Element::GllIntegrationWeightForOrder(options.PolynomialOrder());
    mIntegrationWeightsEta = Element::GllIntegrationWeightForOrder(options.PolynomialOrder());
    mClosureMapping = Element::ClosureMapping(options.PolynomialOrder(), NumberDimensions());

    mNumberIntegrationPointsEps = mIntegrationCoordinatesEps.size();
    mNumberIntegrationPointsEta = mIntegrationCoordinatesEta.size();
    mNumberIntegrationPoints = mNumberIntegrationPointsEps * mNumberIntegrationPointsEta;

    // Allocate element matrices.
    mElementForce.setZero(mNumberIntegrationPoints);
    mElementStrain.setZero(2, mNumberIntegrationPoints);
    mElementDisplacement.setZero(mNumberIntegrationPoints);
    mIntegratedStiffnessMatrix.setZero(mNumberIntegrationPoints);

}


void Acoustic::registerFieldVectors(Mesh *mesh) {

    mesh->registerFieldVector("displacement", true, false);
    mesh->registerFieldVector("acceleration", false, false);
    mesh->registerFieldVector("velocity", false, false);
    mesh->registerFieldVector("force", false, true);

}

void Acoustic::constructStiffnessMatrix(Mesh *mesh) {

    // Test
    for (int i = 0; i < mNumberIntegrationPoints; i++) { mElementDisplacement(i) = mIntegrationCoordinatesEta[i%5]; }

    mesh->getFieldOnElement(mElementDisplacement, mClosureMapping, "displacement", LocalElementNumber());

    int itr = 0;
    Eigen::VectorXd divergence(mNumberIntegrationPoints);
    Eigen::VectorXd external_forcing(mNumberIntegrationPoints);
    Eigen::Matrix<double,2,1> test_function_gradient;
    Eigen::Matrix<double,2,2> inverse_Jacobian;
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            // Eps and eta coordinates.
            double eta = mIntegrationCoordinatesEta[eta_index];
            double eps = mIntegrationCoordinatesEps[eps_index];

            // Get and invert Jacobian.
            inverse_Jacobian = jacobianAtPoint(eps, eta).inverse();

            // Calculate strain. Save for kernel calculations.
            mElementStrain(0,itr) = mGradientOperator.row(eps_index).dot(
                    epsVectorStride(mElementDisplacement, eta_index));
            mElementStrain(1,itr) = mGradientOperator.row(eta_index).dot(
                    etaVectorStride(mElementDisplacement, eta_index));
            mElementStrain.col(itr) = inverse_Jacobian * mElementStrain.col(itr);

            // Calculate test function derivatives.
            test_function_gradient(0) = mGradientOperator.row(eps_index).sum();
            test_function_gradient(1) = mGradientOperator.row(eta_index).sum();
            test_function_gradient = inverse_Jacobian * test_function_gradient;

            // Get material parameters at this node.
            double velocity = interpolateShapeFunctions(eps, eta).dot(mMaterialVelocityAtVertices);
            divergence(itr) = velocity * test_function_gradient.dot(mElementStrain.col(itr));
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
                                         (Square::mIntegrationWeightsEps(eps_index) *
                                                 mIntegrationWeightsEta(eta_index) * determinant) *
                        source->fire(mTime);
            }

            itr++;
        }
    }

    mElementForce = external_forcing - mIntegratedStiffnessMatrix;
    mesh->setFieldOnElement(mElementForce, mClosureMapping, "force", LocalElementNumber());

}

void Acoustic::interpolateMaterialProperties(ExodusModel &model) {

    mMaterialVelocityAtVertices = __interpolateMaterialProperties(model, "velocity");
    mMaterialDensityAtVertices = __interpolateMaterialProperties(model, "density");

}

double Acoustic::evaluateShapeFunctions(const double &eps, const double &eta, const int &itr) {

    Eigen::VectorXd test(mNumberIntegrationPoints);
    interpolate_order4_square(eps, eta, test.data());

    return test(itr);

}