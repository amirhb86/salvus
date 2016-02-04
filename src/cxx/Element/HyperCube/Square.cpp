//
// Created by Michael Afanasiev on 2016-01-30.
//

#include <iostream>
#include <petscdm.h>
#include <petscdmplex.h>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include "Square.h"

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
    for (int i = 0; i < mNumberVertex; i++) {
        mVertexCoordinates(0,i) = coordinates_element[mapping_to_reference_element[mNumberDimensions*i+0]];
        mVertexCoordinates(1,i) = coordinates_element[mapping_to_reference_element[mNumberDimensions*i+1]];
    }

}

void Square::attachIntegrationPoints() {

    assert(NumberIntegrationPoints() == mNumberIntegrationPointsEps*mNumberIntegrationPointsEta);

    int point = 0;
    std::vector<PetscReal> ni(mVertexCoordinates.size());
    mIntegrationPoints.resize(mNumberIntegrationPoints * mNumberDimensions);
    for (auto eta: mIntegrationCoordinatesEta) {
        for (auto eps: mIntegrationCoordinatesEps) {

            Eigen::Vector4d coefficients = interpolateShapeFunctions(eps, eta);

            mIntegrationPoints[point + 0] += coefficients.dot(mVertexCoordinates.row(0));
            mIntegrationPoints[point + 1] += coefficients.dot(mVertexCoordinates.row(1));
            point += mNumberDimensions;

        }
    }
}

PetscReal Square::n0(const PetscReal eps, const PetscReal eta) {
    return 0.25 * (1.0 - eps) * (1.0 - eta);
}

PetscReal Square::n1(const PetscReal eps, const PetscReal eta) {
    return 0.25 * (1.0 + eps) * (1.0 - eta);
}

PetscReal Square::n2(const PetscReal eps, const PetscReal eta) {
    return 0.25 * (1.0 - eps) * (1.0 + eta);
}

PetscReal Square::n3(const PetscReal eps, const PetscReal eta) {
    return 0.25 * (1.0 + eps) * (1.0 + eta);
}

PetscReal Square::dn0deps(const PetscReal eta) {
    return (-1) * (1 - eta) / 4.0;
}

PetscReal Square::dn1deps(const PetscReal eta) {
    return (+1) * (1 - eta) / 4.0;
}

PetscReal Square::dn2deps(const PetscReal eta) {
    return (-1) * (1 + eta) / 4.0;
}

PetscReal Square::dn3deps(const PetscReal eta) {
    return (+1) * (1 + eta) / 4.0;
}

PetscReal Square::dn0deta(const PetscReal eps) {
    return (1 - eps) * -1.0 / 4.0;
}

PetscReal Square::dn1deta(const PetscReal eps) {
    return (1 + eps) * -1.0 / 4.0;
}

PetscReal Square::dn2deta(const PetscReal eps) {
    return (1 - eps) * 1.0 / 4.0;
}

PetscReal Square::dn3deta(const PetscReal eps) {
    return (1 + eps) * 1.0 / 4.0;

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

    Eigen::Vector4d material_at_vertices(mNumberVertex);
    for (auto i = 0; i < mNumberVertex; i++) {
        material_at_vertices(i) = model.getMaterialParameterAtPoint({mVertexCoordinates(0, i),
                                                                     mVertexCoordinates(1, i)},
                                                                    parameter_name);
    }
    return material_at_vertices;

}

void Square::attachSource(std::vector<Source*> sources) {

    for (auto &source: sources) {
        if (mCheckHull(source->LocationX(), source->LocationZ())) {
            std::cout << "HEY!" << std::endl;
            mSources.push_back(source);
        }
    }

}

bool Square::mCheckHull(double x, double z) {

    int n_neg = 0;
    int n_pos = 0;
    Eigen::Vector2d test_point; test_point << x, z;
    for (auto i = 0; i < mNumberVertex; i++) {
        Eigen::Vector2d p0 = mVertexCoordinates.col((i + 0) % mNumberVertex);
        Eigen::Vector2d p1 = mVertexCoordinates.col((i + 1) % mNumberVertex);
        std::cout << p1 << std::endl;
        Eigen::Vector2d v_seg = p1 - p0;
        Eigen::Vector2d p_seg = test_point - p0;
        double x_0 = v_seg(0) * p_seg(1) - v_seg(1) * p_seg(0);
        if (x_0 <= 0) n_neg++;
        if (x_0 >= 0) n_pos++;
    }
    std::cout << n_neg << ' ' << n_pos << std::endl;

    return n_neg == mNumberVertex || n_pos == mNumberVertex ? true : false;

}
