//
// Created by Michael Afanasiev on 2016-01-30.
//

#include "Acoustic.h"

Acoustic::Acoustic(Options options): Square(options) {

    // Allocate element vectors.
    mElementStress.setZero(mNumberIntegrationPoints);
    mIntegratedSource.setZero(mNumberIntegrationPoints);
    mElementDisplacement.setZero(mNumberIntegrationPoints);
    mIntegratedStiffnessMatrix.setZero(mNumberIntegrationPoints);

    // Strain matrix.
    mElementStrain.setZero(2, mNumberIntegrationPoints);

}

void Acoustic::computeStiffnessTerm() {

    // Test
    for (int i = 0; i < mNumberIntegrationPoints; i++) { mElementDisplacement(i) = mIntegrationCoordinatesEta[i%5]; }

    int itr = 0;
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
            mElementStress(itr) = velocity * test_function_gradient.dot(mElementStrain.col(itr));
            itr++;

        }
    }

    itr = 0;
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            double eps = mIntegrationCoordinatesEps[eps_index];
            double eta = mIntegrationCoordinatesEps[eta_index];
            mIntegratedStiffnessMatrix(itr) = jacobianAtPoint(eps, eta).determinant() *
                    (Square::mIntegrationWeightsEps.dot(epsVectorStride(mElementStress, eps_index)) +
                     mIntegrationWeightsEta.dot(etaVectorStride(mElementStress, eps_index)));

            itr++;
        }
    }

}

void Acoustic::interpolateMaterialProperties(ExodusModel &model) {

    mMaterialVelocityAtVertices = __interpolateMaterialProperties(model, "velocity");
    mMaterialDensityAtVertices = __interpolateMaterialProperties(model, "density");

}

void Acoustic::computeSourceTerm() {

    if (! mSources.size()) { mIntegratedSource.setZero(); return; }

    for (auto &source: mSources) {
        interpolate_order4_square(source->ReferenceLocationEps(), source->ReferenceLocationEta(),
                                  mIntegratedSource.data());
        for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
            for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {
                double eps = mIntegrationCoordinatesEps[eps_index];
                double eta = mIntegrationCoordinatesEta[eta_index];
                mIntegratedSource[eps_index + eta_index*mNumberIntegrationPointsEps] /=
                        (mIntegrationWeightsEps(eps_index) * mIntegrationWeightsEta(eta_index) *
                        jacobianAtPoint(eps, eta).determinant()) * source->fire(mTime);
            }
        }
    }
}

void Acoustic::checkOutFields(Mesh *mesh) {

    mElementDisplacement = mesh->getFieldOnElement(
            (int) AcousticFields::displacement , mElementNumber, mClosureMapping);

}

void Acoustic::checkInField(Mesh *mesh) {

    mesh->setFieldOnElement((int) AcousticFields::force, mElementNumber, mClosureMapping, mElementDisplacement);

}

void Acoustic::computeSurfaceTerm() {



}
