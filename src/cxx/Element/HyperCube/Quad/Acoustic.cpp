//
// Created by Michael Afanasiev on 2016-01-30.
//

#include "Acoustic.h"
#include "Autogen/order4_square.h"

Acoustic::Acoustic(Options options): Quad(options) {

    // Allocate element vectors.
    mMassMatrix.setZero(mNumberIntegrationPoints);
    mIntegratedSource.setZero(mNumberIntegrationPoints);
    mElementDisplacement.setZero(mNumberIntegrationPoints);
    mIntegratedStiffnessMatrix.setZero(mNumberIntegrationPoints);

    // Strain matrix.
    mElementStrain.setZero(2, mNumberIntegrationPoints);
    mElementStress.setZero(2, mNumberIntegrationPoints);

}

void Acoustic::computeStiffnessTerm() {

    // Test
//    int j = 0;
//    for (auto i = 0; i < mElementDisplacement.size(); i++) { if (!(i%5) && i>0) j++; mElementDisplacement[i] = mIntegrationCoordinatesEps[j]; }

    int itr = 0;
    Eigen::VectorXd jacobian_determinant(mNumberIntegrationPoints);
    Eigen::Vector2d epsStrain;
    Eigen::Vector2d etaStrain;
    Eigen::Matrix<double,2,1> test_function_gradient;
    Eigen::Matrix<double,2,2> inverse_Jacobian;
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            // Eps and eta coordinates.
            double eta = mIntegrationCoordinatesEta[eta_index];
            double eps = mIntegrationCoordinatesEps[eps_index];

            // Get and invert Jacobian.
            jacobian_determinant(itr) = jacobianAtPoint(eps, eta).determinant();
            inverse_Jacobian = jacobianAtPoint(eps, eta).inverse();

            // Calculate strain. Save for kernel calculations.
//            epsStrain()
            mElementStrain(0,itr) = mGradientOperator.row(eps_index).dot(
                    epsVectorStride(mElementDisplacement, eta_index));
            mElementStrain(1,itr) = mGradientOperator.row(eta_index).dot(
                    etaVectorStride(mElementDisplacement, eta_index));
            mElementStrain.col(itr) = inverse_Jacobian * mElementStrain.col(itr);

            // Get material parameters at this node.
            double velocity = interpolateShapeFunctions(eps, eta).dot(mMaterialVelocityAtVertices);
            mElementStress.col(itr) = mElementStrain.col(itr) * velocity;
            itr++;

        }
    }

    itr = 0;
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            double eps = mIntegrationCoordinatesEps[eps_index];
            double eta = mIntegrationCoordinatesEps[eta_index];

            Eigen::VectorXd jac_stride_eps = epsVectorStride(jacobian_determinant, eta_index);
            Eigen::VectorXd str_stride_eps = epsVectorStride(mElementStress.row(0), eta_index);
            Eigen::VectorXd deriv_eps = mGradientOperator.col(eps_index);

            Eigen::VectorXd jac_stride_eta = etaVectorStride(jacobian_determinant, eta_index);
            Eigen::VectorXd str_stride_eta = etaVectorStride(mElementStress.row(1), eta_index);
            Eigen::VectorXd deriv_eta = mGradientOperator.col(eta_index);

            mIntegratedStiffnessMatrix(itr) = mIntegrationWeightsEta[eta_index] * mIntegrationWeightsEps.dot(
                            (jac_stride_eps.array() * str_stride_eps.array() * deriv_eps.array()).matrix()) +
                    mIntegrationWeightsEps[eps_index] * mIntegrationWeightsEta.dot(
                            (jac_stride_eta.array() * str_stride_eta.array() * deriv_eta.array()).matrix());

            itr++;
        }
    }


}

void Acoustic::interpolateMaterialProperties(ExodusModel &model) {

    mMaterialVelocityAtVertices = __interpolateMaterialProperties(model, "velocity");
    mMaterialDensityAtVertices = __interpolateMaterialProperties(model, "density");

}

void Acoustic::computeSourceTerm() {

    // Check that this integrates to the value of the source (delta function).

    mIntegratedSource.setZero();
    if (! mSources.size()) { return; }

    Eigen::VectorXd current_source(mNumberIntegrationPoints);
    for (auto &source: mSources) {
        interpolate_order4_square(source->ReferenceLocationEps(), source->ReferenceLocationEta(),
                                  current_source.data());
        for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
            for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {
                double eps = mIntegrationCoordinatesEps[eps_index];
                double eta = mIntegrationCoordinatesEta[eta_index];
                current_source[eps_index + eta_index*mNumberIntegrationPointsEps] /=
                        (mIntegrationWeightsEps(eps_index) * mIntegrationWeightsEta(eta_index)) *
                        jacobianAtPoint(eps, eta).determinant();
            }
        }

        current_source *= source->fire(mTime);
        mIntegratedSource(12) = source->fire(mTime);//current_source;
        std::cout << "SOURCE " << mIntegratedSource.maxCoeff() << std::endl;
    }

}

void Acoustic::checkOutFields(Mesh *mesh) {

    mElementDisplacement = mesh->getFieldOnElement(
            (int) AcousticFields::displacement , mElementNumber, mClosureMapping);

}

void Acoustic::checkInField(Mesh *mesh) {

    mesh->setFieldOnElement((int) AcousticFields::force, mElementNumber, mClosureMapping,
                            mIntegratedSource - mIntegratedStiffnessMatrix);

}

void Acoustic::computeSurfaceTerm() {

    if (MPI::COMM_WORLD.Get_rank()) return;
    Eigen::VectorXd test(25);
    interpolate_order4_square(-0.4, -0.7, test.data());
    interpolate_eta_derivative_order4_square(0.7, 0.4, test.data());

}

void Acoustic::assembleMassMatrix() {

    Eigen::VectorXd tmp(mNumberIntegrationPoints);
    double density = mMaterialDensityAtVertices.mean();
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            double eps = mIntegrationCoordinatesEps[eps_index];
            double eta = mIntegrationCoordinatesEta[eta_index];
            diagonal_mass_matrix_order4_square(eps, eta, density, tmp.data());
            mMassMatrix += tmp * mIntegrationWeightsEps[eps_index] * mIntegrationWeightsEta[eta_index] *
                    jacobianAtPoint(eps, eta).determinant();

        }
    }

}
