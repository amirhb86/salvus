//
// Created by Michael Afanasiev on 2016-01-30.
//

#include "Acoustic.h"

AcousticQuad::AcousticQuad(Options options): Quad(options) {

    // Allocate element vectors.
    mMassMatrix.setZero(mNumberIntegrationPoints);

    // Strain matrix.
    mElementStrain.setZero(2, mNumberIntegrationPoints);

}

Eigen::MatrixXd AcousticQuad::computeStiffnessTerm(const Eigen::MatrixXd &displacement) {

    // Current gll point.
    int itr = 0;

    // Data structures we'll need here. Static arrays are allocated for free.
    // TODO: Look into a better way to deal with the temporarily allocated vectors. I believe that due to
    // TODO: RVO the integratedStiffnessMatrix should be ok, but perhaps jacobian_determinant could be handled better.
    Eigen::Matrix<double,2,2> inverse_Jacobian;
    Eigen::MatrixXd velocity_gradient(2, mNumberIntegrationPoints);
    Eigen::VectorXd detJ(mNumberIntegrationPoints);
    Eigen::VectorXd integratedStiffnessMatrix(mNumberIntegrationPoints);
    Eigen::Matrix<double,2,2> Jinv;
    double detJi;

    // Loop over all GLL points once to calculate the gradient of the pressure (u).
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

            // Eps and eta coordinates.
            double eta = mIntegrationCoordinatesEta[eta_index];
            double eps = mIntegrationCoordinatesEps[eps_index];

            // Get Jacobian determinant and its inverse
            std::tie(Jinv,detJi) = inverseJacobianAtPoint(eps,eta);
            detJ[itr] = detJi;

            // Calculate gradient. Save for kernel calculations.
            mElementStrain(0,itr) = mGradientOperator.row(eps_index).dot(
                    epsVectorStride(displacement, eta_index));
            mElementStrain(1,itr) = mGradientOperator.row(eta_index).dot(
                    etaVectorStride(displacement, eps_index));

            // Get material parameters at this node.
            double velocity = interpolateShapeFunctions(eps, eta).dot(mMaterialVelocityAtVertices);

            // apply material velocity (note v^2)
            velocity_gradient.col(itr) = mElementStrain.col(itr) * velocity * velocity;

            // map gradient from reference element to actual element coordinates.
            velocity_gradient.col(itr) = Jinv * velocity_gradient.col(itr);

            itr++;

        }
    }

    // Loop over all gll points again. Apply shape function derivative and integrate.
    itr = 0;
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {
            double eta = mIntegrationCoordinatesEta[eta_index];
            double eps = mIntegrationCoordinatesEps[eps_index];

            Eigen::Matrix<double,2,2> Jinv;
            double detJi;
            std::tie(Jinv,detJi) = inverseJacobianAtPoint(eps,eta);


            // map reference gradient (lr,ls) to this element (lx,lz)
            auto lr = mGradientOperator.col(eps_index);
            auto ls = mGradientOperator.col(eta_index);
            Eigen::VectorXd lx(mNumberIntegrationPointsEps);
            Eigen::VectorXd lz(mNumberIntegrationPointsEta);
            for(int i=0;i<mNumberIntegrationPointsEta;i++) {
                Eigen::VectorXd lrsi(2);
                lrsi[0] = lr[i];
                lrsi[1] = ls[i];
                auto lxzi = Jinv*lrsi;
                lx[i] = lxzi[0];
                lz[i] = lxzi[1];
            }

            integratedStiffnessMatrix(itr) =
                mIntegrationWeightsEta(eta_index) *
                mIntegrationWeightsEps.dot(((epsVectorStride(detJ, eta_index)).array() *
                                            (epsVectorStride(velocity_gradient.row(0), eta_index)).array() *
                                            (lx).array()).matrix()) +
                mIntegrationWeightsEps(eps_index) *
                mIntegrationWeightsEta.dot(((etaVectorStride(detJ, eps_index)).array() *
                                            (etaVectorStride(velocity_gradient.row(1), eps_index)).array() *
                                            (lz).array()).matrix());
            
            itr++;
            
            // old version
            // EQ. (35) in Komatitsch, 2002.
            // TODO: What do you think of this? Too concise?
            // integratedStiffness_rs(0) = mIntegrationWeightsEta(eta_index) *
            //     mIntegrationWeightsEps.dot((
            //                                 (epsVectorStride(jacobian_determinant, eta_index)).array() *
            //                                 (epsVectorStride(velocity_gradient.row(0), eta_index)).array() *
            //                                 (mGradientOperator.col(eps_index)).array()).matrix());

            // integratedStiffness_rs(1) = mIntegrationWeightsEps(eps_index) *
            //     mIntegrationWeightsEta.dot((
            //                                 (etaVectorStride(jacobian_determinant, eps_index)).array() *
            //                                 (etaVectorStride(velocity_gradient.row(1), eps_index)).array() *
            //                                 (mGradientOperator.col(eta_index)).array()).matrix());
            
            
            // integratedStiffnessMatrix(itr) = integratedStiffness_rs(0) + integratedStiffness_rs(1);
            
            

            
        }
    }

    return integratedStiffnessMatrix;
}

void AcousticQuad:: interpolateMaterialProperties(ExodusModel *model) {

    // Vp (m/s).
    mMaterialVelocityAtVertices = __interpolateMaterialProperties(model, "VP");

}

Eigen::MatrixXd AcousticQuad::computeSourceTerm(double time) {

    // Initialize source vector (note: due to RVO I believe no memory re-allocation is occuring).
    Eigen::VectorXd F = Eigen::VectorXd::Zero(mNumberIntegrationPoints);

    // For all sources tagging along with this element.
    for (auto &source: mSources) {

        // TODO: May make this more efficient (i.e. allocation every loop)?
        // Evaluate shape functions at source (eps, eta). Save the lagrange coefficients in current_source.
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
                    (mIntegrationWeightsEps(eps_index) * mIntegrationWeightsEta(eta_index) * detJ);


            }
        }

        // Scale by the source amplitude.
        current_source *= source->fire(time);
        F += current_source;

    }

    return F;
}

void AcousticQuad::computeSurfaceTerm() {

    std::cout << mElementNumber << std::endl;

}

void AcousticQuad::assembleElementMassMatrix(Mesh *mesh) {

    int i=0;
    Eigen::Matrix<double,2,2> Jinv;
    double detJ;
    Eigen::VectorXd elementMassMatrix(mNumberIntegrationPoints);
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {
            double eps = mIntegrationCoordinatesEps[eps_index];
            double eta = mIntegrationCoordinatesEta[eta_index];
            std::tie(Jinv,detJ) = inverseJacobianAtPoint(eps,eta);
            elementMassMatrix[i] = detJ * mIntegrationWeightsEps[eps_index] * mIntegrationWeightsEta[eta_index];
            i++;
        }
    }
    // assemble to shared nodes
    mesh->addFieldFromElement("m", mElementNumber, mClosureMapping, elementMassMatrix);
}

void AcousticQuad::setInitialCondition(Mesh* mesh, Eigen::VectorXd& pts_x,Eigen::VectorXd& pts_z) {
    
    double Lx = 10.0;
    double Lz = 10.0;
    const double PI = std::atan(1.0)*4.0;
    Eigen::VectorXd un = (PI/Lx*(pts_x.array()-Lx/2.0)).sin()*(PI/Lz*(pts_z.array()-Lz/2)).sin();
    Eigen::VectorXd vn = 0*pts_x;
    Eigen::VectorXd an = 0*pts_x;    
    mesh->setFieldFromElement("displacement", mElementNumber, mClosureMapping, un);
    mesh->setFieldFromElement("velocity", mElementNumber, mClosureMapping, vn);
    mesh->setFieldFromElement("acceleration_", mElementNumber, mClosureMapping, an);
    
}

Eigen::VectorXd AcousticQuad::exactSolution(Eigen::VectorXd& pts_x,Eigen::VectorXd& pts_z,double time) {

    const double PI = std::atan(1.0)*4.0;
    double Lx = 10.0;
    double Lz = 10.0;
    Eigen::VectorXd un_xz = (PI/Lx*(pts_x.array()-Lx/2.0)).sin()*(PI/Lz*(pts_z.array()-Lz/2)).sin();
    double velocity = mMaterialVelocityAtVertices.mean();
    double un_t = cos(PI/Lx*sqrt(2)*time*velocity);
    return un_t*un_xz;
}
