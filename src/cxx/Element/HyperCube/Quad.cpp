//
// Created by Michael Afanasiev on 2016-01-30.
//

#include <mpi.h>
#include <iostream>
#include <petscdm.h>
#include <petscdmplex.h>
#include <tuple>
#include "Quad.h"
#include "Quad/Autogen/quad_autogen.h"
#include "Quad/Acoustic.h"
#include "Quad/Elastic.h"

/*
 * STATIC variables WHICH ARE ONLY ON THE REFERENCE ELEMENT.
 */
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

// Eigen::VectorXi Quad::FaceClosureMapping(const int order, const int dimension) {
//     Eigen::VectorXi face_closure_mapping(order+1);
//     if (dimension == 2) {
//         if (order == 4) {
//             face_closure_mapping <<
//                 3, 0, 1, 2, 4;
//         }
//     }
//     return face_closure_mapping;
// }


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
    
    for (int i = 0; i < mNumberVertex; i++) {        
        mVertexCoordinates(0,i) = coordinates_element[mNumberDimensions*i+0];
        mVertexCoordinates(1,i) = coordinates_element[mNumberDimensions*i+1];
    }

    // Save element center
    mElementCenter << mVertexCoordinates.row(0).mean(),
            mVertexCoordinates.row(1).mean();

}

std::tuple<Eigen::Matrix2d, PetscReal> Quad::inverseJacobianAtPoint(PetscReal eps, PetscReal eta) {

    double v1x = mVertexCoordinates(0,0);
    double v2x = mVertexCoordinates(0,1);
    double v3x = mVertexCoordinates(0,2);
    double v4x = mVertexCoordinates(0,3);
    double v1z = mVertexCoordinates(1,0);
    double v2z = mVertexCoordinates(1,1);
    double v3z = mVertexCoordinates(1,2);
    double v4z = mVertexCoordinates(1,3);
    
    double r = (eps+1.0);
    double s = (eta+1.0);
    double detJ = -r*v1x*v3z/8 + r*v1x*v4z/8 + r*v1z*v3x/8 - r*v1z*v4x/8 + r*v2x*v3z/8 -
            r*v2x*v4z/8 - r*v2z*v3x/8 + r*v2z*v4x/8 - s*v1x*v2z/8 + s*v1x*v3z/8 +
            s*v1z*v2x/8 - s*v1z*v3x/8 - s*v2x*v4z/8 + s*v2z*v4x/8 + s*v3x*v4z/8 -
            s*v3z*v4x/8 + v1x*v2z/4 - v1x*v4z/4 - v1z*v2x/4 + v1z*v4x/4 + v2x*v4z/4 -
            v2z*v4x/4;
    
    // J        =   [dx/dr, dy/dr;
    //               dx/ds, dy/ds]
    // J^{-1}   =   [dr/dx ds/dx;
    //               dr/dz ds/dz]
    double rx = (1/detJ)*(r*(v1z - v2z)/4 + r*(v3z - v4z)/4 - v1z/2 + v4z/2);
    double rz = (1/detJ)*(-r*(v1x - v2x)/4 - r*(v3x - v4x)/4 + v1x/2 - v4x/2);
    double sx = (1/detJ)*(-s*(v1z - v2z + v3z - v4z)/4 + v1z/2 - v2z/2);
    double sz = (1/detJ)*(s*(v1x - v2x + v3x - v4x)/4 - v1x/2 + v2x/2);

    Eigen::Matrix<double,2,2> inverseJacobian;
    inverseJacobian << rx, sx,
                       rz, sz;
    
    return std::make_tuple(inverseJacobian, detJ);
}

Eigen::Vector4d Quad::__interpolateMaterialProperties(ExodusModel *model, std::string parameter_name) {

    Eigen::Vector4d material_at_vertices(mNumberVertex);

    for (auto i = 0; i < mNumberVertex; i++) {
        material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(
                mElementCenter, parameter_name, i);
    }
    return material_at_vertices;

}

void Quad::attachSource(std::vector<Source*> sources) {

    for (auto &source: sources) {
        if (mCheckHull(source->PhysicalLocationX(), source->PhysicalLocationZ())) {
            Eigen::Vector2d reference_location = inverseCoordinateTransform(source->PhysicalLocationX(),
                                                                            source->PhysicalLocationZ());
            source->setReferenceLocationEps(reference_location(0));
            source->setReferenceLocationEta(reference_location(1));
            mSources.push_back(source);
        }
    }

}

bool Quad::mCheckHull(double x, double z) {
    int n_neg = 0;
    int n_pos = 0;
    std::vector<int> edge_mapping {0, 1, 2, 3, 0};
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

Eigen::Vector2d Quad::inverseCoordinateTransform(const double &x_real, const double &z_real) {

    double v1x = mVertexCoordinates.row(0)[0];
    double v2x = mVertexCoordinates.row(0)[1];
    double v3x = mVertexCoordinates.row(0)[2];
    double v4x = mVertexCoordinates.row(0)[3];
    double v1z = mVertexCoordinates.row(1)[0];
    double v2z = mVertexCoordinates.row(1)[1];
    double v3z = mVertexCoordinates.row(1)[2];
    double v4z = mVertexCoordinates.row(1)[3];
    
    // Using Newton iterations
    // https://en.wikipedia.org/wiki/Newton%27s_method#Nonlinear_systems_of_equations
    // J_F(xn)(x_{n+1} - x_n) = -F(x_n)
    // Solve for x_{n+1}
    // where J_F(x_n) is jacobian:
    // https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
    double tol = 1e-6;
    int num_iter = 0;
    Eigen::Vector2d solution {0.0, 0.0}; // initial guess at (0.0,0.0)
    while (true) {

        double r = solution(0);
        double s = solution(1);

        Eigen::Matrix2d jacobian;
        // mapping from reference quad [-1,1]x[-1,1] to *this* element
        double Tx = v1x+(v2x-v1x)*(r+1)/2 + (v4x+(v3x-v4x)*(r+1)/2 - v1x - (v2x-v1x)*(r+1)/2)*(s+1)/2;
        double Tz = v1z+(v2z-v1z)*(r+1)/2 + (v4z+(v3z-v4z)*(r+1)/2 - v1z - (v2z-v1z)*(r+1)/2)*(s+1)/2;
        Eigen::Vector2d objective_function {x_real - Tx, z_real - Tz};

        // see element_matrices.py
        // J = [-v1x/2 + v2x/2 + (s + 1)*(v1x/2 - v2x/2 + v3x/2 - v4x/2)/2, -v1x/2 + v4x/2 - (r + 1)*(-v1x + v2x)/4 + (r + 1)*(v3x - v4x)/4],
        // [-v1z/2 + v2z/2 + (s + 1)*(v1z/2 - v2z/2 + v3z/2 - v4z/2)/2, -v1z/2 + v4z/2 - (r + 1)*(-v1z + v2z)/4 + (r + 1)*(v3z - v4z)/4]])
        
        jacobian << (-v1x/2 + v2x/2 + (s + 1)*(v1x/2 - v2x/2 + v3x/2 - v4x/2)/2), -v1x/2 + v4x/2 - (r + 1)*(-v1x + v2x)/4 + (r + 1)*(v3x - v4x)/4, -v1z/2 + v2z/2 + (s + 1)*(v1z/2 - v2z/2 + v3z/2 - v4z/2)/2, -v1z/2 + v4z/2 - (r + 1)*(-v1z + v2z)/4 + (r + 1)*(v3z - v4z)/4;
                    
        if ((objective_function.array().abs() < tol).all()) {
            return solution;
        } else {
            solution += (jacobian.inverse() * objective_function);
        }
        if(num_iter > 10) {
            std::cerr << "inverseCoordinateTransform: TOO MANY ITERATIONS!\n";
            exit(1);
        }
        num_iter++;
        
    }
    
    

}

// TODO: Maybe should be moved to the constructor and made static.
void Quad::setupGradientOperator() {

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
    mNumberDimensions = 2;
    mNumberVertex = 4;
    mPolynomialOrder = options.PolynomialOrder();

    // mVertexCoordinates has 4 vertices
    mVertexCoordinates.resize(2,mNumberVertex);

    // Gll points.
    mNumberDofVertex = 1;
    mNumberDofEdge = mPolynomialOrder - 1;
    mNumberDofFace = (mPolynomialOrder - 1) * (mPolynomialOrder - 1);

    // Integration points.
    mIntegrationCoordinatesEps = Quad::GllPointsForOrder(options.PolynomialOrder());
    mIntegrationCoordinatesEta = Quad::GllPointsForOrder(options.PolynomialOrder());
    mIntegrationWeightsEps = Quad::GllIntegrationWeightForOrder(options.PolynomialOrder());
    mIntegrationWeightsEta = Quad::GllIntegrationWeightForOrder(options.PolynomialOrder());
    mClosureMapping = Quad::ClosureMapping(options.PolynomialOrder(), mNumberDimensions);
    // mFaceClosureMapping = Quad::FaceClosureMapping(options.PolynomialOrder(), mNumberDimensions);

    // Save number of integration points.
    mNumberIntegrationPointsEps = mIntegrationCoordinatesEps.size();
    mNumberIntegrationPointsEta = mIntegrationCoordinatesEta.size();
    mNumberIntegrationPoints = mNumberIntegrationPointsEps * mNumberIntegrationPointsEta;

    // setup evaluated derivatives of test functions
    setupGradientOperator();
}

// global x-z points on all nodes
std::tuple<Eigen::VectorXd,Eigen::VectorXd> Quad::buildNodalPoints() {
		
  assert(mNumberIntegrationPoints == mNumberIntegrationPointsEps*mNumberIntegrationPointsEta);
	
  std::vector<PetscReal> ni(mVertexCoordinates.size());
	
  Eigen::VectorXd nodalPoints_x(mNumberIntegrationPoints);
  Eigen::VectorXd nodalPoints_z(mNumberIntegrationPoints);

  int idx=0;
  for(auto i = 0; i < mNumberIntegrationPointsEta; i++) {	
      for(auto j = 0; j < mNumberIntegrationPointsEps; j++) {
	
          double eps = mIntegrationCoordinatesEps(j);
          double eta = mIntegrationCoordinatesEta(i);

          // reference mapping below uses [0,1]x[0,1] reference square
          double r = (eps+1)/2;
          double s = (eta+1)/2;

          // assumes right-hand rule vertex layout (same as PETSc)
          double v1x = mVertexCoordinates.row(0)[0];
          double v2x = mVertexCoordinates.row(0)[1];
          double v3x = mVertexCoordinates.row(0)[2];
          double v4x = mVertexCoordinates.row(0)[3];
          double v1z = mVertexCoordinates.row(1)[0];
          double v2z = mVertexCoordinates.row(1)[1];
          double v3z = mVertexCoordinates.row(1)[2];
          double v4z = mVertexCoordinates.row(1)[3];

          nodalPoints_x(idx) = v1x+(v2x-v1x)*r + (v4x+(v3x-v4x)*r - v1x - (v2x-v1x)*r)*s;
          nodalPoints_z(idx) = v1z+(v2z-v1z)*r + (v4z+(v3z-v4z)*r - v1z - (v2z-v1z)*r)*s;

          idx++;
      }
  }

  return std::make_tuple(nodalPoints_x,nodalPoints_z);
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

double Quad::integrateField(const Eigen::VectorXd &field) {

    double val = 0;
    Eigen::Matrix<double,2,2> inverse_Jacobian;
    double detJ;
    for (int i = 0; i < mNumberIntegrationPointsEta; i++) {
        for (int j = 0; j < mNumberIntegrationPointsEps; j++) {

            double eps = mIntegrationCoordinatesEps(j);
            double eta = mIntegrationCoordinatesEta(i);
            std::tie(inverse_Jacobian,detJ) = inverseJacobianAtPoint(eps,eta);
            val += field(j+i*mNumberIntegrationPointsEps) * mIntegrationWeightsEps(j) *
                mIntegrationWeightsEta(i) * detJ;

        }
    }

    return val;

}
