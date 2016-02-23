//
// Created by Michael Afanasiev on 2016-01-30.
//

#include "Acoustic.h"

Acoustic::Acoustic(Options options): Quad(options) {

    // Allocate element vectors.
    mMassMatrix.setZero(mNumberIntegrationPoints);

    // Strain matrix.
    mElementStrain.setZero(2, mNumberIntegrationPoints);

}

Eigen::MatrixXd Acoustic::computeStiffnessTerm(const Eigen::MatrixXd &displacement) {

    // Current gll point.
    int itr = 0;

    // Data structures we'll need here. Static arrays are allocated for free.
    // TODO: Look into a better way to deal with the temporarily allocated vectors. I believe that due to
    // TODO: RVO the integratedStiffnessMatrix should be ok, but perhaps jacobian_determinant could be handled better.
    Eigen::Matrix<double,2,2> inverse_Jacobian;
    Eigen::VectorXd stress(2, mNumberIntegrationPoints);
    Eigen::VectorXd jacobian_determinant(mNumberIntegrationPoints);
    Eigen::VectorXd integratedStiffnessMatrix(mNumberIntegrationPoints);

    // Loop over all GLL points once to calculate the stress.
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            // Eps and eta coordinates.
            double eta = mIntegrationCoordinatesEta[eta_index];
            double eps = mIntegrationCoordinatesEps[eps_index];

            // Get and invert Jacobian.
            jacobian_determinant(itr) = jacobianAtPoint(eps, eta).determinant();
            inverse_Jacobian = jacobianAtPoint(eps, eta).inverse();

            // Calculate strain. Save for kernel calculations.
            mElementStrain(0,itr) = mGradientOperator.row(eps_index).dot(
                    epsVectorStride(displacement, eta_index));
            mElementStrain(1,itr) = mGradientOperator.row(eta_index).dot(
                    etaVectorStride(displacement, eps_index));
            mElementStrain.col(itr) = inverse_Jacobian * mElementStrain.col(itr);

            // Get material parameters at this node.
            double velocity = interpolateShapeFunctions(eps, eta).dot(mMaterialVelocityAtVertices);

            // Stress = material parameter * strain.
            stress.col(itr) = mElementStrain.col(itr) * velocity * velocity;

            // Convert stress to physical co-ordinates.
            stress.col(itr) = inverse_Jacobian * stress.col(itr);

            itr++;

        }
    }

    // Loop over all gll points again. Apply shape function derivative and integrate.
    itr = 0;
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            // Eq. (35) in Komatitsch, 2002.
            // TODO: What do you think of this? Too concise?
            integratedStiffnessMatrix(itr) =
                    mIntegrationWeightsEta(eta_index) *
                    mIntegrationWeightsEps.dot((
                            (epsVectorStride(jacobian_determinant, eta_index)).array() *
                            (epsVectorStride(stress.row(0), eta_index)).array() *
                            (mGradientOperator.col(eps_index)).array()).matrix()) +

                    mIntegrationWeightsEps(eps_index) *
                    mIntegrationWeightsEta.dot((
                            (etaVectorStride(jacobian_determinant, eps_index)).array() *
                            (etaVectorStride(stress.row(1), eps_index)).array() *
                            (mGradientOperator.col(eta_index)).array()).matrix());

            itr++;
        }
    }

    return integratedStiffnessMatrix;
}

void Acoustic::interpolateMaterialProperties(ExodusModel *model) {

    mMaterialVelocityAtVertices = __interpolateMaterialProperties(model, "velocity");
    mMaterialDensityAtVertices = __interpolateMaterialProperties(model, "density");

}

Eigen::MatrixXd Acoustic::computeSourceTerm() {

    // Initialize source vector (note: due to RVO I believe no memory re-allocation is occuring).
    Eigen::VectorXd F = Eigen::VectorXd::Zero(mNumberIntegrationPoints);

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
        F(12) += source->fire(mTime);

    }

    return F;
}

void Acoustic::computeSurfaceTerm() {


}

void Acoustic::assembleElementMassMatrix(Mesh *mesh) {

    int i=0;
    Eigen::VectorXd elementMassMatrix(mNumberIntegrationPoints);
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {
            double eps = mIntegrationCoordinatesEps[eps_index];
            double eta = mIntegrationCoordinatesEta[eta_index];
            elementMassMatrix[i] = mIntegrationWeightsEps[eps_index] * mIntegrationWeightsEta[eta_index] *
                    jacobianAtPoint(eps, eta).determinant();
            i++;
        }
    }
    // assemble to shared nodes
    mesh->addFieldFromElement("mass_matrix", mElementNumber, mClosureMapping, elementMassMatrix);
}

void Acoustic::setInitialCondition(Mesh* mesh, Eigen::VectorXd& pts_x,Eigen::VectorXd& pts_z) {
    
    double Lx = 10.0;
    double Lz = 10.0;
    const double PI = std::atan(1.0)*4.0;
//    Eigen::VectorXd un = (PI/Lx*(pts_x.array()-Lx/2.0)).sin()*(PI/Lz*(pts_z.array()-Lz/2)).sin();
    Eigen::VectorXd un = (PI/Lx*(pts_x.array()-Lx/2.0)).sin();//*(PI/Lz*(pts_z.array()-Lz/2)).sin();
    mesh->setFieldFromElement("displacement", mElementNumber, mClosureMapping, un);
}

