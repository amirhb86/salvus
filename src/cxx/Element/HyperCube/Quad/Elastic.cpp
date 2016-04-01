//
// Created by Michael Afanasiev on 2016-02-17.
//

#include "Elastic.h"

// Elemental fields definition.
const std::vector<std::string> mElementalFields {"ux", "uy"};

Elastic::Elastic(Options options): Quad(options) {

    mMassMatrix.setZero(mNumberIntegrationPoints);
    mElementStrain.setZero(2, mNumberIntegrationPoints*2);

}

Eigen::MatrixXd Elastic::computeSourceTerm(double time) {

    // Initialize source vector (note: due to RVO I believe no memory re-allocation is occuring).
    Eigen::MatrixXd F = Eigen::MatrixXd::Zero(mNumberIntegrationPoints, 2);

    // For all sources tagging along with this element.
    for (auto &source: mSources) {

        // TODO: May make this more efficient (i.e. allocation every loop)
        Eigen::VectorXd current_source = interpolateLagrangePolynomials(
                source->ReferenceLocationEps(), source->ReferenceLocationEta(), mPolynomialOrder);

        // Loop over gll points
        for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
            for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

                double eps = mIntegrationCoordinatesEps[eps_index];
                double eta = mIntegrationCoordinatesEta[eta_index];

                double detJ;
                Eigen::Matrix2d Jinv;
                std::tie(Jinv, detJ) = inverseJacobianAtPoint(eps,eta);

                // Calculate the coefficients needed to integrate to the delta function.
                current_source[eps_index + eta_index*mNumberIntegrationPointsEps] /=
                        (mIntegrationWeightsEps(eps_index) * mIntegrationWeightsEta(eta_index) *
                        detJ);

            }
        }

        // Scale by the source amplitude.
        current_source *= source->fire(time);

        // TODO: Right now this is hardcoded for a source in the x-direction.
        F.col(0) = F.col(0) + current_source;
    }

    return F;

}

void Elastic::computeSurfaceTerm() {

}

void Elastic::assembleElementMassMatrix(Mesh *mesh) {

    int i = 0;
    Eigen::VectorXd elementMassMatrix(mNumberIntegrationPoints);
    double density = mRhoAtVertices.mean();
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            Eigen::Matrix2d jInv;
            double jDet;
            double eps = mIntegrationCoordinatesEps[eps_index];
            double eta = mIntegrationCoordinatesEta[eta_index];
            std::tie(jInv, jDet) = inverseJacobianAtPoint(eps,eta);
            elementMassMatrix[i] = density * mIntegrationWeightsEps[eps_index] * mIntegrationWeightsEta[eta_index] *
                                   jDet;
            i++;

        }
    }

    mesh->addFieldFromElement("m", mElementNumber, mClosureMapping, elementMassMatrix);

}

Eigen::MatrixXd Elastic::computeStiffnessTerm(const Eigen::MatrixXd &displacement) {

    int itr = 0;
    Eigen::Matrix2d inverse_jacobian, temp_stress;
    Eigen::VectorXd jacobian_determinant(mNumberIntegrationPoints);
    Eigen::VectorXd element_stress_xx(mNumberIntegrationPoints);
    Eigen::VectorXd element_stress_xz(mNumberIntegrationPoints);
    Eigen::VectorXd element_stress_zx(mNumberIntegrationPoints);
    Eigen::VectorXd element_stress_zz(mNumberIntegrationPoints);
    Eigen::MatrixXd integratedStiffnessMatrix(mNumberIntegrationPoints,2);
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            // Eps and eta coordinates.
            double eta = mIntegrationCoordinatesEta[eta_index];
            double eps = mIntegrationCoordinatesEps[eps_index];

            std::tie(inverse_jacobian, jacobian_determinant(itr)) = inverseJacobianAtPoint(eps, eta);

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
                    epsVectorStride(displacement.col(0), eta_index));
            mElementStrain(rx,cz) = mGradientOperator.row(eta_index).dot(
                    etaVectorStride(displacement.col(0), eps_index));
            mElementStrain(rz,cx) = mGradientOperator.row(eps_index).dot(
                    epsVectorStride(displacement.col(1), eta_index));
            mElementStrain(rz,cz) = mGradientOperator.row(eta_index).dot(
                    etaVectorStride(displacement.col(1), eps_index));

            // Convert derivatives to physical co-ordinates through the jacobian.
            mElementStrain.block<2,2>(rx,cx) = inverse_jacobian * mElementStrain.block<2,2>(rx,cx);

            // Get material parameters at this integration point.
            double c11 = interpolateAtPoint(eps, eta).dot(mC11AtVertices);
            double c13 = interpolateAtPoint(eps, eta).dot(mC13AtVertices);
            double c15 = interpolateAtPoint(eps, eta).dot(mC15AtVertices);
            double c33 = interpolateAtPoint(eps, eta).dot(mC33AtVertices);
            double c35 = interpolateAtPoint(eps, eta).dot(mC35AtVertices);
            double c55 = interpolateAtPoint(eps, eta).dot(mC55AtVertices);

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
            integratedStiffnessMatrix(itr,0) =
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

            integratedStiffnessMatrix(itr,1) =
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

void Elastic::attachMaterialProperties(ExodusModel *model) {

    // TODO: C15, C35.
    mRhoAtVertices = __attachMaterialProperties(model, "RHO");
    mC11AtVertices = __attachMaterialProperties(model, "C11");
    mC13AtVertices = __attachMaterialProperties(model, "C13");
//    mC15AtVertices = __attachMaterialProperties(model, "C15");
    mC33AtVertices = __attachMaterialProperties(model, "C33");
//    mC35AtVertices = __attachMaterialProperties(model, "C35");
    mC55AtVertices = __attachMaterialProperties(model, "C55");

    mC31AtVertices = mC13AtVertices = Eigen::MatrixXd::Zero(mNumberVertex,1);
    mC51AtVertices = mC15AtVertices = Eigen::MatrixXd::Zero(mNumberVertex,1);
    mC53AtVertices = mC35AtVertices = Eigen::MatrixXd::Zero(mNumberVertex,1);

}
