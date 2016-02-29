//
// Created by Michael Afanasiev on 2016-01-30.
//

#include <iostream>
#include <petscdm.h>
#include <petscdmplex.h>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <tuple>
#include "Quad.h"
#include "Quad/Autogen/quad_autogen.h"
#include "Quad/Acoustic.h"
#include "Quad/Elastic.h"

/*
 * STATIC FUNCTIONS WHICH ARE ONLY ON THE REFERENCE ELEMENT.
 */
int Quad::mNumberDofVertex;
int Quad::mNumberDofEdge;
int Quad::mNumberDofFace;
int Quad::mNumberDofVolume;
int Quad::mNumberIntegrationPointsEps;
int Quad::mNumberIntegrationPointsEta;
int Quad::mNumberIntegrationPoints;
int Quad::mPolynomialOrder;

Eigen::VectorXi Quad::mClosureMapping;
Eigen::MatrixXd Quad::mGradientOperator;
Eigen::VectorXd Quad::mIntegrationWeightsEps;
Eigen::VectorXd Quad::mIntegrationWeightsEta;
Eigen::VectorXd Quad::mIntegrationCoordinatesEps;
Eigen::VectorXd Quad::mIntegrationCoordinatesEta;

double Quad::n0(const double &eps, const double &eta) { return 0.25 * (1.0 - eps) * (1.0 - eta); }
double Quad::n1(const double &eps, const double &eta) { return 0.25 * (1.0 + eps) * (1.0 - eta); }
double Quad::n2(const double &eps, const double &eta) { return 0.25 * (1.0 - eps) * (1.0 + eta); }
double Quad::n3(const double &eps, const double &eta) { return 0.25 * (1.0 + eps) * (1.0 + eta); }
double Quad::dn0deps(const double &eta) { return (-1) * (1 - eta) / 4.0; }
double Quad::dn1deps(const double &eta) { return (+1) * (1 - eta) / 4.0; }
double Quad::dn2deps(const double &eta) { return (-1) * (1 + eta) / 4.0; }
double Quad::dn3deps(const double &eta) { return (+1) * (1 + eta) / 4.0; }
double Quad::dn0deta(const double &eps) { return (1 - eps) * -1.0 / 4.0; }
double Quad::dn1deta(const double &eps) { return (1 + eps) * -1.0 / 4.0; }
double Quad::dn2deta(const double &eps) { return (1 - eps) * 1.0 / 4.0; }
double Quad::dn3deta(const double &eps) { return (1 + eps) * 1.0 / 4.0; }


Eigen::VectorXd Quad::GllPointsForOrder(const int order) {
    Eigen::VectorXd gll_points(order+1);
    if (order == 1) {
        gll_coordinates_order1_square(gll_points.data());
    } else if (order == 2) {
        gll_coordinates_order2_square(gll_points.data());
    } else if (order == 3) {
        gll_coordinates_order3_square(gll_points.data());
    } else if (order == 4) {
        gll_coordinates_order4_square(gll_points.data());
    } else if (order == 5) {
        gll_coordinates_order5_square(gll_points.data());
    } else if (order == 6) {
        gll_coordinates_order6_square(gll_points.data());
    } else if (order == 7) {
        gll_coordinates_order7_square(gll_points.data());
    } else if (order == 8) {
        gll_coordinates_order8_square(gll_points.data());
    } else if (order == 9) {
        gll_coordinates_order9_square(gll_points.data());
    } else if (order == 10) {
        gll_coordinates_order10_square(gll_points.data());
    }
    return gll_points;
}

Eigen::VectorXd Quad::GllIntegrationWeightForOrder(const int order) {
    Eigen::VectorXd integration_weights(order+1);
    if (order == 1) {
        gll_weights_order1_square(integration_weights.data());
    } else if (order == 2) {
        gll_weights_order2_square(integration_weights.data());
    } else if (order == 3) {
        gll_weights_order3_square(integration_weights.data());
    } else if (order == 4) {
        gll_weights_order4_square(integration_weights.data());
    } else if (order == 5) {
        gll_weights_order5_square(integration_weights.data());
    } else if (order == 6) {
        gll_weights_order6_square(integration_weights.data());
    } else if (order == 7) {
        gll_weights_order7_square(integration_weights.data());
    } else if (order == 8) {
        gll_weights_order8_square(integration_weights.data());
    } else if (order == 9) {
        gll_weights_order9_square(integration_weights.data());
    } else if (order == 10) {
        gll_weights_order10_square(integration_weights.data());
    }
    return integration_weights;
}

Eigen::VectorXi Quad::ClosureMapping(const int order, const int dimension) {
    Eigen::VectorXi closure_mapping((order+1)*(order+1));
    if (order == 1) {
        closure_mapping_order1_square(closure_mapping.data());
    } else if (order == 2) {
        closure_mapping_order2_square(closure_mapping.data());
    } else if (order == 3) {
        closure_mapping_order3_square(closure_mapping.data());
    } else if (order == 4) {
        closure_mapping_order4_square(closure_mapping.data());
    } else if (order == 5) {
        closure_mapping_order5_square(closure_mapping.data());
    } else if (order == 6) {
        closure_mapping_order6_square(closure_mapping.data());
    } else if (order == 7) {
        closure_mapping_order7_square(closure_mapping.data());
    } else if (order == 8) {
        closure_mapping_order8_square(closure_mapping.data());
    } else if (order == 9) {
        closure_mapping_order9_square(closure_mapping.data());
    } else if (order == 10) {
        closure_mapping_order10_square(closure_mapping.data());
    }
    return closure_mapping;
}


Eigen::Map<const Eigen::VectorXd> Quad::epsVectorStride(const Eigen::VectorXd &function,
                                                        const int &eta_index) {
    return Eigen::Map<const Eigen::VectorXd> (
            function.data() + eta_index * mNumberIntegrationPointsEta,
            mNumberIntegrationPointsEps);
}

Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>> Quad::etaVectorStride(const Eigen::VectorXd &function,
                                                                                 const int &eta_index) {
    return Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>> (
            function.data() + eta_index, mNumberIntegrationPointsEta,
            Eigen::InnerStride<> (mNumberIntegrationPointsEps));
}


Eigen::Vector4d Quad::interpolateShapeFunctions(const double &eps, const double &eta) {
    Eigen::Vector4d coefficients;
    coefficients(0) = n0(eps, eta);
    coefficients(1) = n1(eps, eta);
    coefficients(2) = n2(eps, eta);
    coefficients(3) = n3(eps, eta);
    return coefficients;
}

void Quad::attachVertexCoordinates(DM &distributed_mesh) {

    Vec coordinates_local;
    PetscInt coordinate_buffer_size;
    PetscSection coordinate_section;
    PetscReal *coordinates_buffer = NULL;

    DMGetCoordinatesLocal(distributed_mesh, &coordinates_local);
    DMGetCoordinateSection(distributed_mesh, &coordinate_section);
    DMPlexVecGetClosure(distributed_mesh, coordinate_section, coordinates_local, mElementNumber,
                        &coordinate_buffer_size, &coordinates_buffer);
    std::vector<PetscReal> coordinates_element(coordinates_buffer, coordinates_buffer+coordinate_buffer_size);
    DMPlexVecRestoreClosure(distributed_mesh, coordinate_section, coordinates_local, mElementNumber,
                            &coordinate_buffer_size, &coordinates_buffer);

    // Reorder to desired vertex ordering.
    std::vector<double> vertex_coordinates_ordered ;
    std::vector<PetscInt> mapping_to_reference_element {6, 7, 0, 1, 4, 5, 2, 3};
    for (int i = 0; i < mNumberVertex; i++) {
        mVertexCoordinates(0,i) = coordinates_element[mapping_to_reference_element[mNumberDimensions*i+0]];
        mVertexCoordinates(1,i) = coordinates_element[mapping_to_reference_element[mNumberDimensions*i+1]];
    }

}

Eigen::Matrix<double,2,2> Quad::jacobianAtPoint(PetscReal eps, PetscReal eta) {
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

Eigen::Vector4d Quad::__interpolateMaterialProperties(ExodusModel *model, std::string parameter_name) {

    Eigen::Vector4d material_at_vertices(mNumberVertex);

    // THIS IS STUPID JUST BECAUSE DENSITY IS NOT IN THE EXODUS FILE YET.
    if (parameter_name == "density") {
        material_at_vertices.setConstant(3.0);
        return material_at_vertices;
    }

    for (auto i = 0; i < mNumberVertex; i++) {
        material_at_vertices(i) = model->getMaterialParameterAtPoint({mVertexCoordinates(0, i),
                                                                     mVertexCoordinates(1, i)},
                                                                    parameter_name);
    }
    return material_at_vertices;

}

void Quad::attachSource(std::vector<Source*> sources) {

    for (auto &source: sources) {
        if (mCheckHull(source->PhysicalLocationX(), source->PhysicalLocationZ())) {
            Eigen::Vector2d reference_location = inverseCoordinateTransform(source->PhysicalLocationX(),
                                                                            source->PhysicalLocationZ(),
                                                                            0.0, 0.0);
            source->setReferenceLocationEps(reference_location(0));
            source->setReferenceLocationEta(reference_location(1));
            mSources.push_back(source);
        }
    }

}

bool Quad::mCheckHull(double x, double z) {
    int n_neg = 0;
    int n_pos = 0;
    std::vector<int> edge_mapping {0, 1, 3, 2, 0};
    Eigen::Vector2d test_point; test_point << x, z;
    for (auto i = 0; i < mNumberVertex; i++) {
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
    return n_neg == mNumberVertex || n_pos == mNumberVertex;
}

Eigen::Vector2d Quad::inverseCoordinateTransform(const double &x_real, const double &z_real,
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

// TODO: Maybe should be moved to the constructor and made static.
void Quad::readGradientOperator() {

    double eta = mIntegrationCoordinatesEta[0];
    mGradientOperator.resize(mNumberIntegrationPointsEta, mNumberIntegrationPointsEps);
    Eigen::MatrixXd test(mNumberIntegrationPointsEta, mNumberIntegrationPointsEps);
    for (auto i=0; i < mNumberIntegrationPointsEps; i++) {
        double eps = mIntegrationCoordinatesEps[i];
        if (mPolynomialOrder == 1) {
            interpolate_eps_derivative_order1_square(eta, test.data());
        } else if (mPolynomialOrder == 2) {
            interpolate_eps_derivative_order2_square(eps, eta, test.data());
        } else if (mPolynomialOrder == 3) {
            interpolate_eps_derivative_order3_square(eps, eta, test.data());
        } else if (mPolynomialOrder == 4) {
            interpolate_eps_derivative_order4_square(eps, eta, test.data());
        } else if (mPolynomialOrder == 5) {
            interpolate_eps_derivative_order5_square(eps, eta, test.data());
        } else if (mPolynomialOrder == 6) {
            interpolate_eps_derivative_order6_square(eps, eta, test.data());
        } else if (mPolynomialOrder == 7) {
            interpolate_eps_derivative_order7_square(eps, eta, test.data());
        } else if (mPolynomialOrder == 8) {
            interpolate_eps_derivative_order8_square(eps, eta, test.data());
        } else if (mPolynomialOrder == 9) {
            interpolate_eps_derivative_order9_square(eps, eta, test.data());
        } else if (mPolynomialOrder == 10) {
            interpolate_eps_derivative_order10_square(eps, eta, test.data());
        }
        mGradientOperator.row(i) = test.col(0);
    }

}

Quad::Quad(Options options) {

    // Basic properties.
    mPolynomialOrder = options.PolynomialOrder();

    // Gll points.
    mNumberDofVolume = 0;
    mNumberDofVertex = 1;
    mNumberDofEdge = mPolynomialOrder - 1;
    mNumberDofFace = (mPolynomialOrder - 1) * (mPolynomialOrder - 1);

    // Integration points.
    mIntegrationCoordinatesEps = Quad::GllPointsForOrder(options.PolynomialOrder());
    mIntegrationCoordinatesEta = Quad::GllPointsForOrder(options.PolynomialOrder());
    mIntegrationWeightsEps = Quad::GllIntegrationWeightForOrder(options.PolynomialOrder());
    mIntegrationWeightsEta = Quad::GllIntegrationWeightForOrder(options.PolynomialOrder());
    mClosureMapping = Quad::ClosureMapping(options.PolynomialOrder(), mNumberDimensions);

    // Save number of integration points.
    mNumberIntegrationPointsEps = mIntegrationCoordinatesEps.size();
    mNumberIntegrationPointsEta = mIntegrationCoordinatesEta.size();
    mNumberIntegrationPoints = mNumberIntegrationPointsEps * mNumberIntegrationPointsEta;

}

// complete assembly between common nodes
// void Quad::assembleMassMatrix(Mesh *mesh) {
//     mesh->setFieldOnElement("mass_matrix", mElementNumber, mClosureMapping, mMassMatrix);
// }

// global x-z points on all nodes
std::tuple<Eigen::VectorXd,Eigen::VectorXd> Quad::buildNodalPoints(Mesh* mesh) {
		
  assert(mNumberIntegrationPoints == mNumberIntegrationPointsEps*mNumberIntegrationPointsEta);
	
  std::vector<PetscReal> ni(mVertexCoordinates.size());
	
  Eigen::VectorXd nodalPoints_x(mNumberIntegrationPoints);
  Eigen::VectorXd nodalPoints_z(mNumberIntegrationPoints);

  int idx=0;
  for(auto i = 0; i < mNumberIntegrationPointsEta; i++) {	
      for(auto j = 0; j < mNumberIntegrationPointsEps; j++) {
	
          double eps = mIntegrationCoordinatesEps(j);
          double eta = mIntegrationCoordinatesEta(i);
		
          nodalPoints_x(idx) += interpolateShapeFunctions(eps, eta).dot(mVertexCoordinates.row(0));      
          nodalPoints_z(idx) += interpolateShapeFunctions(eps, eta).dot(mVertexCoordinates.row(1));
          idx++;
	
      }
  }
  // push nodal locations to shared dofs
  mesh->setFieldFromElement("nodes_x", mElementNumber, mClosureMapping, nodalPoints_x);
  mesh->setFieldFromElement("nodes_z", mElementNumber, mClosureMapping, nodalPoints_z);

  return std::make_tuple(nodalPoints_x,nodalPoints_z);
    
}

Eigen::VectorXd Quad::checkOutFieldElement(Mesh *mesh, const std::string name) {

    return mesh->getFieldOnElement(name, mElementNumber, mClosureMapping);

}

void Quad::checkInFieldElement(Mesh *mesh, const Eigen::VectorXd &field, const std::string name) {

    mesh->setFieldOnElement(name, mElementNumber, mClosureMapping, field);

}

Quad *Quad::factory(Options options) {

    std::string physics(options.PhysicsSystem());
    try {
        if (physics == "acoustic") {
            return new Acoustic(options);
        } else if (physics == "elastic") {
            return new Elastic(options);
        } else {
            throw std::runtime_error("Runtime Error: Element physics " + physics + " not supported");
        }
    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI::COMM_WORLD.Abort(-1);
        return nullptr;
    }
}

Eigen::VectorXd Quad::interpolateLagrangePolynomials(const double eps, const double eta, const int p_order) {

    assert(p_order > 0 && p_order < 11);

    int n_points = (p_order + 1) * (p_order + 1);
    Eigen::VectorXd gll_coeffs(n_points);
    if (p_order == 1) {
        interpolate_order1_square(eps, eta, gll_coeffs.data());
    } else if (p_order == 2) {
        interpolate_order2_square(eps, eta, gll_coeffs.data());
    } else if (p_order == 3) {
        interpolate_order3_square(eps, eta, gll_coeffs.data());
    } else if (p_order == 4) {
        interpolate_order4_square(eps, eta, gll_coeffs.data());
    } else if (p_order == 5) {
        interpolate_order5_square(eps, eta, gll_coeffs.data());
    } else if (p_order == 6) {
        interpolate_order6_square(eps, eta, gll_coeffs.data());
    } else if (p_order == 7) {
        interpolate_order7_square(eps, eta, gll_coeffs.data());
    } else if (p_order == 8) {
        interpolate_order8_square(eps, eta, gll_coeffs.data());
    } else if (p_order == 9) {
        interpolate_order9_square(eps, eta, gll_coeffs.data());
    } else if (p_order == 10) {
        interpolate_order10_square(eps, eta, gll_coeffs.data());
    }
    return gll_coeffs;
}
