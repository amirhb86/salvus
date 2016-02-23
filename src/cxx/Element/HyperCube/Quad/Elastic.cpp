//
// Created by Michael Afanasiev on 2016-02-17.
//

#include "Elastic.h"

Elastic::Elastic(Options options): Quad(options) {

    mMassMatrix.setZero(mNumberIntegrationPoints);
    mElementStrain.setZero(2, mNumberIntegrationPoints*2);

}

void Elastic::checkInField(Mesh *mesh) {

}

Eigen::MatrixXd Elastic::computeSourceTerm() {

    // Initialize source vector (note: due to RVO I believe no memory re-allocation is occuring).
    Eigen::MatrixXd F = Eigen::MatrixXd::Zero(2, mNumberIntegrationPoints);

    // For all sources tagging along with this element.
    for (auto &source: mSources) {

        // TODO: May make this more efficient (i.e. allocation every loop)
        Eigen::VectorXd current_source(mNumberIntegrationPoints);

        // Evaluate shape functions at source (eps, eta). Save the lagrange coefficients in current_source.
        interpolate_order4_square(source->ReferenceLocationEps(), source->ReferenceLocationEta(),
                                  current_source.data());

        // Loop over gll points
        for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
            for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

                double eps = mIntegrationCoordinatesEps[eps_index];
                double eta = mIntegrationCoordinatesEta[eta_index];

                // Calculate the coefficients needed to integrate to the delta function.
                current_source[eps_index + eta_index*mNumberIntegrationPointsEps] /=
                        (mIntegrationWeightsEps(eps_index) * mIntegrationWeightsEta(eta_index)) *
                        jacobianAtPoint(eps, eta).determinant();

            }
        }

        // Scale by the source amplitude.
        current_source *= source->fire(mTime);

        // TODO: Add current source to F. Right now, this isn't working quite right, so I'm just putting it at
        // a gll point.
        F(0,12) += source->fire(mTime);
    }

    return F;

}

void Elastic::computeSurfaceTerm() {

}

void Elastic::assembleElementMassMatrix(Mesh *mesh) {

    int i = 0;
    Eigen::VectorXd elementMassMatrix(mNumberIntegrationPoints);
    double density = 3.0;//mRhoAtVertices.mean();
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            double eps = mIntegrationCoordinatesEps[eps_index];
            double eta = mIntegrationCoordinatesEta[eta_index];
            elementMassMatrix[i] = density * mIntegrationWeightsEps[eps_index] * mIntegrationWeightsEta[eta_index] *
                                   jacobianAtPoint(eps, eta).determinant();
            i++;

        }
    }

    mesh->addFieldFromElement("mass_matrix", mElementNumber, mClosureMapping, elementMassMatrix);

}

Eigen::MatrixXd Elastic::computeStiffnessTerm(const Eigen::MatrixXd &displacement) {

    int itr = 0;
    Eigen::Matrix2d inverse_jacobian, temp_stress;
    Eigen::VectorXd jacobian_determinant(mNumberIntegrationPoints);
    Eigen::VectorXd element_stress_xx(mNumberIntegrationPoints);
    Eigen::VectorXd element_stress_xz(mNumberIntegrationPoints);
    Eigen::VectorXd element_stress_zx(mNumberIntegrationPoints);
    Eigen::VectorXd element_stress_zz(mNumberIntegrationPoints);
    Eigen::MatrixXd integratedStiffnessMatrix(2,mNumberIntegrationPoints);
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            // Eps and eta coordinates.
            double eta = mIntegrationCoordinatesEta[eta_index];
            double eps = mIntegrationCoordinatesEps[eps_index];

            jacobian_determinant(itr) = jacobianAtPoint(eps, eta).determinant();
            inverse_jacobian = jacobianAtPoint(eps, eta).inverse();

            /*
             * | uxx   uxz  |
             * |            |
             * | uzx   uzz  |
             * Here the strain is calculated for each gll point. At the moment, the full strain matrix is calculated
             * (see above). It is stored in the pre-allocated matrix containing the strain on each element, in 2x2
             * blocks. This is why there's the 2*itr+0,1 below. */
            int rx = 0;
            int rz = 1;
            int cx = 2*itr+0;
            int cz = 2*itr+1;
            mElementStrain(rx,cx) = mGradientOperator.row(eps_index).dot(
                    epsVectorStride(displacement.row(0), eta_index));
            mElementStrain(rx,cz) = mGradientOperator.row(eta_index).dot(
                    etaVectorStride(displacement.row(0), eps_index));
            mElementStrain(rz,cx) = mGradientOperator.row(eps_index).dot(
                    epsVectorStride(displacement.row(1), eta_index));
            mElementStrain(rz,cz) = mGradientOperator.row(eta_index).dot(
                    etaVectorStride(displacement.row(1), eps_index));

            // Convert derivatives to physical co-ordinates through the jacobian.
            mElementStrain.block<2,2>(rx,cx) = inverse_jacobian * mElementStrain.block<2,2>(rx,cx);

            // Get material parameters at this integration point.
            double c11 = interpolateShapeFunctions(eps, eta).dot(mC11AtVertices);
            double c13 = interpolateShapeFunctions(eps, eta).dot(mC13AtVertices);
            double c15 = interpolateShapeFunctions(eps, eta).dot(mC15AtVertices);
            double c33 = interpolateShapeFunctions(eps, eta).dot(mC33AtVertices);
            double c35 = interpolateShapeFunctions(eps, eta).dot(mC35AtVertices);
            double c55 = interpolateShapeFunctions(eps, eta).dot(mC55AtVertices);

            // Calculate element stress by multiplying strain through by elastic tensor.
            // TODO: Probably don't need this temporary object.
            temp_stress(0,0) = c11 * mElementStrain(rx,cx) + c13 * mElementStrain(rz,cz) +
                               c15 * (mElementStrain(rz,cx) + mElementStrain(rx,cz));
            temp_stress(1,1) = c13 * mElementStrain(rx,cx) + c33 * mElementStrain(rz,cz) +
                               c35 * (mElementStrain(rz,cx) + mElementStrain(rx,cz));
            temp_stress(0,1) = temp_stress(1,0) = c15 * mElementStrain(rx,cx) +
                                                  c35 * mElementStrain(rz,cz) +
                                                  c55 * (mElementStrain(rz,cx) + mElementStrain(rx,cz));

            // Transform stress to physical co-ordinates.
            temp_stress = inverse_jacobian * temp_stress;

            // Copy to nice stress array.
            element_stress_xx(itr) = temp_stress(0,0);
            element_stress_zz(itr) = temp_stress(1,1);
            element_stress_xz(itr) = element_stress_zx(itr) = temp_stress(0,1);

            itr++;

        }
    }

    // Loop over all shape functions again. Apply shape function derivatives and integrate.
    // Stiffness matrix is like kx -> row(0), kz -> row(1).
    itr = 0;
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            /*
             * Eq. (35) in Komatitsch, 2002.
             * TODO: So this is obviously pretty concise. The nice thing is that it recreates equation 35 pretty
             * explicitly. It also (should be) optimized by Eigen... although I'd be pretty impressed if that's true.
             * Anyone know how to read assembly?
             */
            integratedStiffnessMatrix(0,itr) =
                    mIntegrationWeightsEta(eta_index) *
                    mIntegrationWeightsEps.dot(
                            ((epsVectorStride(jacobian_determinant, eta_index)).array() *
                            (epsVectorStride(element_stress_xx, eta_index)).array() *
                            (mGradientOperator.col(eps_index)).array()).matrix()) +

                    mIntegrationWeightsEps(eps_index) *
                    mIntegrationWeightsEta.dot(
                            ((etaVectorStride(jacobian_determinant, eps_index)).array() *
                            (etaVectorStride(element_stress_xz, eps_index)).array() *
                            (mGradientOperator.col(eta_index)).array()).matrix());

            integratedStiffnessMatrix(1,itr) =
                    mIntegrationWeightsEta(eta_index) *
                    mIntegrationWeightsEps.dot(
                            ((epsVectorStride(jacobian_determinant, eta_index)).array() *
                            (epsVectorStride(element_stress_xz, eta_index)).array() *
                             (mGradientOperator.col(eps_index)).array()).matrix()) +

                    mIntegrationWeightsEps(eps_index) *
                    mIntegrationWeightsEta.dot(
                            ((etaVectorStride(jacobian_determinant, eps_index)).array() *
                            (etaVectorStride(element_stress_zz, eps_index)).array() *
                            (mGradientOperator.col(eta_index)).array()).matrix());

            itr++;
        }
    }

    return integratedStiffnessMatrix;
}

void Elastic::interpolateMaterialProperties(ExodusModel *model) {

    mC11AtVertices = __interpolateMaterialProperties(model, "c11");
    mC13AtVertices = __interpolateMaterialProperties(model, "c13");
    mC15AtVertices = __interpolateMaterialProperties(model, "c15");
    mC31AtVertices = __interpolateMaterialProperties(model, "c31");
    mC33AtVertices = __interpolateMaterialProperties(model, "c33");
    mC35AtVertices = __interpolateMaterialProperties(model, "c35");
    mC15AtVertices = __interpolateMaterialProperties(model, "c15");
    mC35AtVertices = __interpolateMaterialProperties(model, "c35");
    mC55AtVertices = __interpolateMaterialProperties(model, "c55");
}

void Elastic::setInitialCondition(Mesh *mesh, Eigen::VectorXd &pts_x, Eigen::VectorXd &ptz_z) {

}