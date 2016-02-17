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

}

void Elastic::computeSurfaceTerm() {

}

void Elastic::assembleElementMassMatrix(Mesh *mesh) {

}

Eigen::MatrixXd Elastic::computeStiffnessTerm(const Eigen::MatrixXd &displacement) {

    int itr = 0;
    Eigen::Matrix2d inverse_jacobian;
    Eigen::VectorXd jacobian_determinant(mNumberIntegrationPoints);
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            // Eps and eta coordinates.
            double eta = mIntegrationCoordinatesEta[eta_index];
            double eps = mIntegrationCoordinatesEps[eps_index];

            jacobian_determinant(itr) = jacobianAtPoint(eps, eta).determinant();
            inverse_jacobian = jacobianAtPoint(eps, eta).inverse();

            /*
             * | uxx   uxy  |
             * |            |
             * | uyx   uyy  |
             * Here the strain is calculated for each gll point. At the moment, the full strain matrix is calculated
             * (see above). It is stored in the pre-allocated matrix containing the strain on each element, in 2x2
             * blocks. This is why there's the 2*itr+0,1 below.
             */
            mElementStrain(0,2*itr+0) = mGradientOperator.row(eps_index).dot(
                    epsVectorStride(displacement.row(0), eta_index));
            mElementStrain(0,2*itr+1) = mGradientOperator.row(eta_index).dot(
                    etaVectorStride(displacement.row(0), eta_index));
            mElementStrain(1,2*itr+0) = mGradientOperator.row(eps_index).dot(
                    epsVectorStride(displacement.row(1), eta_index));
            mElementStrain(1,2*itr+1) = mGradientOperator.row(eta_index).dot(
                    etaVectorStride(displacement.row(1), eta_index));
        }
    }

    Eigen::MatrixXd test;
    return test;
}

void Elastic::interpolateMaterialProperties(ExodusModel *model) {

}

void Elastic::setInitialCondition(Mesh *mesh, Eigen::VectorXd &pts_x, Eigen::VectorXd &ptz_z) {

}
