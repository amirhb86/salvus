//
// Created by Michael Afanasiev on 2016-01-30.
//

#include "Acoustic.h"

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

Eigen::VectorXd Acoustic::computeStiffnessTerm() {

    // Test
//    int j = 0;
//    for (auto i = 0; i < mElementDisplacement.size(); i++) { if (!(i%5) && i>0) j++; mElementDisplacement[i] = mIntegrationCoordinatesEps[j]; }

    int itr = 0;
    Eigen::VectorXd jacobian_determinant(mNumberIntegrationPoints);
    Eigen::Vector2d epsStrain;
    Eigen::Vector2d etaStrain;
    Eigen::Matrix<double,2,1> test_function_gradient;
    Eigen::Matrix<double,2,2> inverse_Jacobian;

    Eigen::VectorXd integratedStiffnessMatrix(mNumberIntegrationPoints);
    
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

            integratedStiffnessMatrix(itr) = mIntegrationWeightsEta[eta_index] * mIntegrationWeightsEps.dot(
                            (jac_stride_eps.array() * str_stride_eps.array() * deriv_eps.array()).matrix()) +
                    mIntegrationWeightsEps[eps_index] * mIntegrationWeightsEta.dot(
                            (jac_stride_eta.array() * str_stride_eta.array() * deriv_eta.array()).matrix());

            itr++;
        }
    }
    return integratedStiffnessMatrix;
}

void Acoustic::interpolateMaterialProperties(ExodusModel *model) {

    mMaterialVelocityAtVertices = __interpolateMaterialProperties(model, "velocity");
    mMaterialDensityAtVertices = __interpolateMaterialProperties(model, "density");

}

Eigen::VectorXd Acoustic::computeSourceTerm() {

    // Check that this integrates to the value of the source (delta function).
    
    mIntegratedSource.setZero();
    Eigen::VectorXd F;
    if(mSources.size() > 0) {
        F.resize(mIntegratedSource.size());
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
            // mIntegratedSource(12) = source->fire(mTime);//current_source;
            // TODO: Why 12?
            F(12) = source->fire(mTime);
        
            std::cout << "SOURCE " << mIntegratedSource.maxCoeff() << std::endl;
            
        }
    }
    return F;
}

void Acoustic::checkOutField(Mesh *mesh) {
    mElementDisplacement = mesh->getFieldOnElement("displacement", mElementNumber, mClosureMapping);
}

// void Acoustic::checkInVector()

void Acoustic::checkInField(Mesh *mesh) {
    auto FminusK = mIntegratedSource - mIntegratedStiffnessMatrix;
    mesh->setFieldOnElement("force", mElementNumber, mClosureMapping,
                            FminusK);
}



void Acoustic::checkInFieldElement(Mesh *mesh,Eigen::VectorXd& field) {
    mesh->setFieldOnElement("force", mElementNumber, mClosureMapping,
                            field);
}

void Acoustic::computeSurfaceTerm() {

    if (MPI::COMM_WORLD.Get_rank()) return;
    Eigen::VectorXd test(25);
    interpolate_order4_square(-0.4, -0.7, test.data());
    interpolate_eta_derivative_order4_square(0.7, 0.4, test.data());

}

void Acoustic::assembleElementMassMatrix(Mesh *mesh) {

    Eigen::VectorXd elementMassMatrix(mNumberIntegrationPoints);
    double density = mMaterialDensityAtVertices.mean();
    int i=0;
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
        for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {
            double eps = mIntegrationCoordinatesEps[eps_index];
            double eta = mIntegrationCoordinatesEta[eta_index];
            elementMassMatrix[i] = density * mIntegrationWeightsEps[eps_index] * mIntegrationWeightsEta[eta_index] *
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
    Eigen::VectorXd un = (PI/Lx*(pts_x.array()-Lx/2.0)).sin()*(PI/Lz*(pts_z.array()-Lz/2)).sin();
    mesh->setFieldFromElement("displacement", mElementNumber, mClosureMapping, un);
}

