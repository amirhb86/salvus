//
// Created by Michael Afanasiev on 2016-01-30.
//

#include <iostream>
#include <petscdm.h>
#include <petscdmplex.h>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include "Square.h"
#include "Square/Autogen/order4_square.h"

PetscReal Square::n0(const PetscReal eps, const PetscReal eta) { return 0.25 * (1.0 - eps) * (1.0 - eta); }
PetscReal Square::n1(const PetscReal eps, const PetscReal eta) { return 0.25 * (1.0 + eps) * (1.0 - eta); }
PetscReal Square::n2(const PetscReal eps, const PetscReal eta) { return 0.25 * (1.0 - eps) * (1.0 + eta); }
PetscReal Square::n3(const PetscReal eps, const PetscReal eta) { return 0.25 * (1.0 + eps) * (1.0 + eta); }
PetscReal Square::dn0deps(const PetscReal eta) { return (-1) * (1 - eta) / 4.0; }
PetscReal Square::dn1deps(const PetscReal eta) { return (+1) * (1 - eta) / 4.0; }
PetscReal Square::dn2deps(const PetscReal eta) { return (-1) * (1 + eta) / 4.0; }
PetscReal Square::dn3deps(const PetscReal eta) { return (+1) * (1 + eta) / 4.0; }
PetscReal Square::dn0deta(const PetscReal eps) { return (1 - eps) * -1.0 / 4.0; }
PetscReal Square::dn1deta(const PetscReal eps) { return (1 + eps) * -1.0 / 4.0; }
PetscReal Square::dn2deta(const PetscReal eps) { return (1 - eps) * 1.0 / 4.0; }
PetscReal Square::dn3deta(const PetscReal eps) { return (1 + eps) * 1.0 / 4.0; }

void Square::attachVertexCoordinates() {

    Vec coordinates_local;
    PetscInt coordinate_buffer_size;
    PetscSection coordinate_section;
    PetscReal *coordinates_buffer = NULL;

    DMGetCoordinatesLocal(mDistributedMesh, &coordinates_local);
    DMGetCoordinateSection(mDistributedMesh, &coordinate_section);
    DMPlexVecGetClosure(mDistributedMesh, coordinate_section, coordinates_local, LocalElementNumber(),
                        &coordinate_buffer_size, &coordinates_buffer);
    std::vector<PetscReal> coordinates_element(coordinates_buffer, coordinates_buffer+coordinate_buffer_size);
    DMPlexVecRestoreClosure(mDistributedMesh, coordinate_section, coordinates_local, LocalElementNumber(),
                            &coordinate_buffer_size, &coordinates_buffer);

    // Reorder to desired vertex ordering.
    std::vector<double> vertex_coordinates_ordered ;
    std::vector<PetscInt> mapping_to_reference_element {6, 7, 0, 1, 4, 5, 2, 3};
    for (int i = 0; i < NumberVertex(); i++) {
        mVertexCoordinates(0,i) = coordinates_element[mapping_to_reference_element[mNumberDimensions*i+0]];
        mVertexCoordinates(1,i) = coordinates_element[mapping_to_reference_element[mNumberDimensions*i+1]];
    }

}

void Square::attachIntegrationPoints() {

    assert(NumberIntegrationPoints() == mNumberIntegrationPointsEps*mNumberIntegrationPointsEta);

    int point = 0;
    std::vector<PetscReal> ni(mVertexCoordinates.size());
    mIntegrationPoints.resize(mNumberIntegrationPoints * mNumberDimensions);
    for (auto i = 0; i < mNumberIntegrationPointsEta; i++) {
        for (auto j = 0; j < mNumberIntegrationPointsEps; j++) {

            double eps = mIntegrationCoordinatesEps(j);
            double eta = mIntegrationCoordinatesEta(i);

            mIntegrationPoints[point + 0] += interpolateShapeFunctions(eps, eta).dot(mVertexCoordinates.row(0));
            mIntegrationPoints[point + 1] += interpolateShapeFunctions(eps, eta).dot(mVertexCoordinates.row(1));
            point += mNumberDimensions;

        }
    }
}

Eigen::Matrix<double,2,2> Square::jacobianAtPoint(PetscReal eps, PetscReal eta) {

    // Set local values.
    Eigen::Matrix<double,2,4> jacobian_multiplier;
    jacobian_multiplier(0,0) = dn0deps(eta);
    jacobian_multiplier(0,1) = dn1deps(eta);
    jacobian_multiplier(0,2) = dn2deps(eta);
    jacobian_multiplier(0,3) = dn3deps(eta);
    jacobian_multiplier(1,0) = dn0deta(eps);
    jacobian_multiplier(1,1) = dn1deta(eps);
    jacobian_multiplier(1,2) = dn2deta(eps);
    jacobian_multiplier(1,3) = dn3deta(eps);

    return jacobian_multiplier * mVertexCoordinates.transpose();

}

Eigen::Vector4d Square::interpolateShapeFunctions(PetscReal eps, PetscReal eta) {

    Eigen::Vector4d coefficients;
    coefficients(0) = n0(eps, eta);
    coefficients(1) = n1(eps, eta);
    coefficients(2) = n2(eps, eta);
    coefficients(3) = n3(eps, eta);
    return coefficients;

}

Eigen::Map<Eigen::VectorXd> Square::epsVectorStride(
        Eigen::VectorXd &function, int &eta_index) {
    return Eigen::Map<Eigen::VectorXd> (

            function.data() + eta_index * mNumberIntegrationPointsEta,
            mNumberIntegrationPointsEps);

}

Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<>> Square::etaVectorStride(
        Eigen::VectorXd &function, int &eta_index) {

    return Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<>> (
            function.data() + eta_index, mNumberIntegrationPointsEta,
            Eigen::InnerStride<> (mNumberIntegrationPointsEps));

}

Eigen::Vector4d Square::__interpolateMaterialProperties(ExodusModel &model, std::string parameter_name) {

    Eigen::Vector4d material_at_vertices(NumberVertex());

    // THIS IS STUPID JUST BECAUSE DENSITY IS NOT IN THE EXODUS FILE YET.
    if (parameter_name == "density") {
        material_at_vertices.setConstant(3.0);
        return material_at_vertices;
    }


    for (auto i = 0; i < NumberVertex(); i++) {
        material_at_vertices(i) = model.getMaterialParameterAtPoint({mVertexCoordinates(0, i),
                                                                     mVertexCoordinates(1, i)},
                                                                    parameter_name);
    }
    return material_at_vertices;

}

void Square::attachSource(std::vector<Source*> sources) {

    for (auto &source: sources) {
        if (mCheckHull(source->PhysicalLocationX(), source->PhysicalLocationZ())) {
            Eigen::Vector2d reference_location = inverseCoordinateTransform(source->PhysicalLocationX(),
                                                                            source->PhysicalLocationZ(),
                                                                            0.0, 0.0);
            source->setReferenceLocationEps(reference_location(0));
            source->setReferenceLocationEta(reference_location(1));
            mSources.push_back(source);
            SetContainsSource(true);
        }
    }

}

bool Square::mCheckHull(double x, double z) {
    int n_neg = 0;
    int n_pos = 0;
    std::vector<int> edge_mapping {0, 1, 3, 2, 0};
    Eigen::Vector2d test_point; test_point << x, z;
    for (auto i = 0; i < NumberVertex(); i++) {
        Eigen::Vector2d p0 = mVertexCoordinates.col(edge_mapping[i+0]);
        Eigen::Vector2d p1 = mVertexCoordinates.col(edge_mapping[i+1]);
        Eigen::Vector2d v_seg = p1 - p0;
        Eigen::Vector2d p_seg = test_point - p0;
        double x_0 = v_seg(0) * p_seg(1) - v_seg(1) * p_seg(0);
        if (x_0 <= 0) {
            n_neg++;
        } else {
            n_pos++;
        }
    }
    return n_neg == NumberVertex() || n_pos == NumberVertex();
}

Eigen::Vector2d Square::inverseCoordinateTransform(const double &x_real, const double &z_real,
                                                   double eps, double eta) {

    double tol = 1e-6;
    Eigen::Vector2d solution {eps, eta};
    while (true) {

        eps = solution(0);
        eta = solution(1);

        Eigen::Matrix2d jacobian;
        Eigen::Vector4d shape_functions {n0(eps, eta), n1(eps, eta), n2(eps, eta), n3(eps, eta)};
        Eigen::Vector4d dNdEps {dn0deps(eta), dn1deps(eta), dn2deps(eta), dn3deps(eta)};
        Eigen::Vector4d dNdEta {dn0deta(eps), dn1deta(eps), dn2deta(eps), dn3deta(eps)};
        Eigen::Vector2d objective_function {x_real - shape_functions.dot(mVertexCoordinates.row(0)),
                                            z_real - shape_functions.dot(mVertexCoordinates.row(1))};

        jacobian << -1 * dNdEps.dot(mVertexCoordinates.row(0)),
                    -1 * dNdEta.dot(mVertexCoordinates.row(0)),
                    -1 * dNdEps.dot(mVertexCoordinates.row(1)),
                    -1 * dNdEta.dot(mVertexCoordinates.row(1));

        if ((objective_function.array().abs() < tol).all()) {
            return solution;
        } else {
            solution -= (jacobian.inverse() * objective_function);
        }

    }

}

void Square::readOperators() {

    int i = 0;
    double eta = mIntegrationCoordinatesEta[0];
    mGradientOperator.resize(mNumberIntegrationPointsEta, mNumberIntegrationPointsEps);
    Eigen::MatrixXd test(mNumberIntegrationPointsEta, mNumberIntegrationPointsEps);
    for (auto i=0; i < mNumberIntegrationPointsEps; i++) {
        double eps = mIntegrationCoordinatesEps[i];
        interpolate_eps_derivative_order4_square(eps, eta, test.data());
        mGradientOperator.row(i) = test.col(0);
        i++;
    }

}
