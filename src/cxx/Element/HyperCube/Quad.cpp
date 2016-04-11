//
// Created by Michael Afanasiev on 2016-01-30.
//

#include <mpi.h>

#include <petscdm.h>
#include <petscdmplex.h>
#include "Quad.h"

double Quad::p1n0(const double &eps, const double &eta) { return 0.25 * (1.0 - eps) * (1.0 - eta); }
double Quad::p1n1(const double &eps, const double &eta) { return 0.25 * (1.0 + eps) * (1.0 - eta); }
double Quad::p1n2(const double &eps, const double &eta) { return 0.25 * (1.0 - eps) * (1.0 + eta); }
double Quad::p1n3(const double &eps, const double &eta) { return 0.25 * (1.0 + eps) * (1.0 + eta); }

Quad::Quad(Options options) {

  // Basic properties.
  mNumDim = 2;
  mNumVtx = 4;
  mPlyOrd = options.PolynomialOrder();

  // mVtxCrd has 4 vertices.
  mElmCtr.resize(2, 1);
  mVtxCrd.resize(2, mNumVtx);

  // Gll points.
  mNumDofVol = 0;
  mNumDofVtx = 1;
  mNumDofEdg = mPlyOrd - 1;
  mNumDofFac = (mPlyOrd - 1) * (mPlyOrd - 1);

  // Integration points.
  mClsMap = Quad::ClosureMapping(options.PolynomialOrder(), mNumDim);
  mIntCrdEps = Quad::GllPointsForOrder(options.PolynomialOrder());
  mIntCrdEta = Quad::GllPointsForOrder(options.PolynomialOrder());
  mIntWgtEps = Quad::GllIntegrationWeightForOrder(options.PolynomialOrder());
  mIntWgtEta = Quad::GllIntegrationWeightForOrder(options.PolynomialOrder());

  // Save number of integration points.
  mNumIntPtsEps = mIntCrdEps.size();
  mNumIntPtsEta = mIntCrdEta.size();
  mNumIntPnt = mNumIntPtsEps * mNumIntPtsEta;

  // setup evaluated derivatives of test functions
  setupGradientOperator();
}

VectorXd Quad::GllPointsForOrder(const int order) {
  VectorXd gll_points(order + 1);
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

VectorXd Quad::GllIntegrationWeightForOrder(const int order) {
  VectorXd integration_weights(order + 1);
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

VectorXi Quad::ClosureMapping(const int order, const int dimension) {
  VectorXi closure_mapping((order + 1) * (order + 1));
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

Map<const VectorXd> Quad::epsVectorStride(const VectorXd &function,
                                          const int &eta_index) {
  return Map<const VectorXd>(
      function.data() + eta_index * mNumIntPtsEta,
      mNumIntPtsEps);
}

Map<const VectorXd, 0, InnerStride<>> Quad::etaVectorStride(const VectorXd &function,
                                                            const int &eta_index) {
  return Map<const VectorXd, 0, InnerStride<>>(
                                               function.data() + eta_index, mNumIntPtsEta,
                                               InnerStride<>(mNumIntPtsEps));
}

void Quad::attachVertexCoordinates(DM &distributed_mesh) {

  Vec coordinates_local;
  PetscInt coordinate_buffer_size;
  PetscSection coordinate_section;
  PetscReal *coordinates_buffer = NULL;

  DMGetCoordinatesLocal(distributed_mesh, &coordinates_local);
  DMGetCoordinateSection(distributed_mesh, &coordinate_section);
  DMPlexVecGetClosure(distributed_mesh, coordinate_section, coordinates_local, mElmNum,
                      &coordinate_buffer_size, &coordinates_buffer);
  std::vector<PetscReal> coordinates_element(coordinates_buffer, coordinates_buffer + coordinate_buffer_size);
  DMPlexVecRestoreClosure(distributed_mesh, coordinate_section, coordinates_local, mElmNum,
                          &coordinate_buffer_size, &coordinates_buffer);

  for (int i = 0; i < mNumVtx; i++) {
    mVtxCrd(0, i) = coordinates_element[mNumDim * i + 0];
    mVtxCrd(1, i) = coordinates_element[mNumDim * i + 1];
  }

  // Save element center
  mElmCtr << mVtxCrd.row(0).mean(),
      mVtxCrd.row(1).mean();

}

std::tuple<Matrix2d, PetscReal> Quad::inverseJacobianAtPoint(
    PetscReal eps, PetscReal eta) {

  double v1x = mVtxCrd(0, 0);
  double v2x = mVtxCrd(0, 1);
  double v3x = mVtxCrd(0, 2);
  double v4x = mVtxCrd(0, 3);
  double v1z = mVtxCrd(1, 0);
  double v2z = mVtxCrd(1, 1);
  double v3z = mVtxCrd(1, 2);
  double v4z = mVtxCrd(1, 3);

  double r = (eps + 1.0);
  double s = (eta + 1.0);
  double detJ = -r * v1x * v3z / 8 + r * v1x * v4z / 8 + r * v1z * v3x / 8 - r * v1z * v4x / 8 + r * v2x * v3z / 8 -
      r * v2x * v4z / 8 - r * v2z * v3x / 8 + r * v2z * v4x / 8 - s * v1x * v2z / 8 + s * v1x * v3z / 8 +
      s * v1z * v2x / 8 - s * v1z * v3x / 8 - s * v2x * v4z / 8 + s * v2z * v4x / 8 + s * v3x * v4z / 8 -
      s * v3z * v4x / 8 + v1x * v2z / 4 - v1x * v4z / 4 - v1z * v2x / 4 + v1z * v4x / 4 + v2x * v4z / 4 -
      v2z * v4x / 4;

  // J        =   [dx/dr, dy/dr;
  //               dx/ds, dy/ds]
  // J^{-1}   =   [dr/dx ds/dx;
  //               dr/dz ds/dz]
  double rx = (1 / detJ) * (r * (v1z - v2z) / 4 + r * (v3z - v4z) / 4 - v1z / 2 + v4z / 2);
  double rz = (1 / detJ) * (-r * (v1x - v2x) / 4 - r * (v3x - v4x) / 4 + v1x / 2 - v4x / 2);
  double sx = (1 / detJ) * (-s * (v1z - v2z + v3z - v4z) / 4 + v1z / 2 - v2z / 2);
  double sz = (1 / detJ) * (s * (v1x - v2x + v3x - v4x) / 4 - v1x / 2 + v2x / 2);

  Matrix2d inverseJacobian;
  inverseJacobian <<  rx, sx,
                      rz, sz;

  return std::make_tuple(inverseJacobian, detJ);
}

Eigen::Vector4d Quad::interpolateAtPoint(double eps, double eta) {
  // see element_metrics.py
  //  interpolation for quadrilaterals based on right-hand rule quad with
  // velocity v0,v1,v2,v3 and eps(x) eta(y,or,z)
  //
  //  v3               vb          v2
  //  +-----------------+----------+ \eps
  //  |                 |          |   ^       \alpha = (eta-(-1))/2
  //  |                 |          |   |       \beta  = (eps-(-1))/2
  //  |              vf + \beta    |   |       va = alpha(v1) + (1-alpha)v0
  //  |                 |  ^       |   |       vb = alpha(v2) + (1-alpha)v3
  //  |                 |  |       |           vf = beta va + (1-beta) vb
  //  |                 |  |       |           group vf by terms v0..v3
  //  |                 |  |       |           = [v0,v1,v2,v3].dot([-0.25*eps*eta - 0.25*eps + 0.25*eta + 0.25,
  //  |              va |  |       |                                0.25*eps*eta + 0.25*eps + 0.25*eta + 0.25,
  //  +-----------------+----------+                                -0.25*eps*eta + 0.25*eps - 0.25*eta + 0.25,
  //  v0 -------------->\alpha    v1  --->\eta                      0.25*eps*eta - 0.25*eps - 0.25*eta + 0.25]
  // ----------------------------------------------------------------------------------------------------------
  Eigen::Vector4d interpolator;
  interpolator <<        
    0.25*eps*eta - 0.25*eps - 0.25*eta + 0.25,
    -0.25*eps*eta + 0.25*eps - 0.25*eta + 0.25,
    0.25*eps*eta + 0.25*eps + 0.25*eta + 0.25,
    -0.25*eps*eta - 0.25*eps + 0.25*eta + 0.25;
  return interpolator;
}

Eigen::Vector4d Quad::__attachMaterialProperties(ExodusModel *model, std::string parameter_name) {

  Eigen::Vector4d material_at_vertices(mNumVtx);

  for (auto i = 0; i < mNumVtx; i++) {
    material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(
        mElmCtr, parameter_name, i);
  }
  return material_at_vertices;

}

void Quad::attachSource(std::vector<Source *> sources) {

  for (auto &source: sources) {
    if (mCheckHull(source->PhysicalLocationX(), source->PhysicalLocationZ())) {
      Vector2d reference_location = inverseCoordinateTransform(source->PhysicalLocationX(),
                                                               source->PhysicalLocationZ());
      source->setReferenceLocationEps(reference_location(0));
      source->setReferenceLocationEta(reference_location(1));
      mSrc.push_back(source);
    }
  }

}

bool Quad::mCheckHull(double x, double z) {
  int n_neg = 0;
  int n_pos = 0;
  std::vector<int> edge_mapping{0, 1, 2, 3, 0};
  Vector2d test_point;
  test_point << x, z;
  for (auto i = 0; i < mNumVtx; i++) {
    Vector2d p0 = mVtxCrd.col(edge_mapping[i + 0]);
    Vector2d p1 = mVtxCrd.col(edge_mapping[i + 1]);
    Vector2d v_seg = p1 - p0;
    Vector2d p_seg = test_point - p0;
    double x_0 = v_seg(0) * p_seg(1) - v_seg(1) * p_seg(0);
    if (x_0 <= 0) {
      n_neg++;
    } else {
      n_pos++;
    }
  }
  return n_neg == mNumVtx || n_pos == mNumVtx;
}

Vector2d Quad::inverseCoordinateTransform(const double &x_real, const double &z_real) {

  double v1x = mVtxCrd.row(0)[0];
  double v2x = mVtxCrd.row(0)[1];
  double v3x = mVtxCrd.row(0)[2];
  double v4x = mVtxCrd.row(0)[3];
  double v1z = mVtxCrd.row(1)[0];
  double v2z = mVtxCrd.row(1)[1];
  double v3z = mVtxCrd.row(1)[2];
  double v4z = mVtxCrd.row(1)[3];

  // Using Newton iterations
  // https://en.wikipedia.org/wiki/Newton%27s_method#Nonlinear_systems_of_equations
  // J_F(xn)(x_{n+1} - x_n) = -F(x_n)
  // Solve for x_{n+1}
  // where J_F(x_n) is jacobian:
  // https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
  double tol = 1e-6;
  int num_iter = 0;
  Vector2d solution{0.0, 0.0}; // initial guess at (0.0,0.0)
  while (true) {

    double r = solution(0);
    double s = solution(1);

    Matrix2d jacobian;

    // mapping from reference quad [-1,1]x[-1,1] to *this* element
    double Tx = v1x + (v2x - v1x) * (r + 1) / 2
        + (v4x + (v3x - v4x) * (r + 1) / 2 - v1x - (v2x - v1x) * (r + 1) / 2) * (s + 1) / 2;
    double Tz = v1z + (v2z - v1z) * (r + 1) / 2
        + (v4z + (v3z - v4z) * (r + 1) / 2 - v1z - (v2z - v1z) * (r + 1) / 2) * (s + 1) / 2;
    Vector2d objective_function{x_real - Tx, z_real - Tz};

    // see element_matrices.py
    // J = [-v1x/2 + v2x/2 + (s + 1)*(v1x/2 - v2x/2 + v3x/2 - v4x/2)/2, -v1x/2 + v4x/2 -
    // (r + 1)*(-v1x + v2x)/4 + (r + 1)*(v3x - v4x)/4],
    // [-v1z/2 + v2z/2 + (s + 1)*(v1z/2 - v2z/2 + v3z/2 - v4z/2)/2,
    // -v1z/2 + v4z/2 - (r + 1)*(-v1z + v2z)/4 + (r + 1)*(v3z - v4z)/4]])

    jacobian << (-v1x / 2 + v2x / 2 + (s + 1) * (v1x / 2 - v2x / 2 + v3x / 2 - v4x / 2) / 2), -v1x / 2 + v4x / 2
        - (r + 1) * (-v1x + v2x) / 4 + (r + 1) * (v3x - v4x) / 4, -v1z / 2 + v2z / 2
        + (s + 1) * (v1z / 2 - v2z / 2 + v3z / 2 - v4z / 2) / 2, -v1z / 2 + v4z / 2 - (r + 1) * (-v1z + v2z) / 4
        + (r + 1) * (v3z - v4z) / 4;

    if ((objective_function.array().abs() < tol).all()) {
      return solution;
    } else {
      solution += (jacobian.inverse() * objective_function);
    }
    if (num_iter > 10) {
      std::cerr << "inverseCoordinateTransform: TOO MANY ITERATIONS!\n";
      exit(1);
    }
    num_iter++;

  }

}

// TODO: Maybe should be moved to the constructor and made static.
void Quad::setupGradientOperator() {

  double eta = mIntCrdEta[0];
  mGrd.resize(mNumIntPtsEta, mNumIntPtsEps);
  Eigen::MatrixXd test(mNumIntPtsEta, mNumIntPtsEps);
  for (auto i = 0; i < mNumIntPtsEps; i++) {
    double eps = mIntCrdEps[i];
    if (mPlyOrd == 1) {
      interpolate_eps_derivative_order1_square(eta, test.data());
    } else if (mPlyOrd == 2) {
      interpolate_eps_derivative_order2_square(eps, eta, test.data());
    } else if (mPlyOrd == 3) {
      interpolate_eps_derivative_order3_square(eps, eta, test.data());
    } else if (mPlyOrd == 4) {
      interpolate_eps_derivative_order4_square(eps, eta, test.data());
    } else if (mPlyOrd == 5) {
      interpolate_eps_derivative_order5_square(eps, eta, test.data());
    } else if (mPlyOrd == 6) {
      interpolate_eps_derivative_order6_square(eps, eta, test.data());
    } else if (mPlyOrd == 7) {
      interpolate_eps_derivative_order7_square(eps, eta, test.data());
    } else if (mPlyOrd == 8) {
      interpolate_eps_derivative_order8_square(eps, eta, test.data());
    } else if (mPlyOrd == 9) {
      interpolate_eps_derivative_order9_square(eps, eta, test.data());
    } else if (mPlyOrd == 10) {
      interpolate_eps_derivative_order10_square(eps, eta, test.data());
    }
    mGrd.row(i) = test.col(0);
  }
}

std::tuple<VectorXd, VectorXd> Quad::buildNodalPoints() {

  assert(mNumIntPnt == mNumIntPtsEps * mNumIntPtsEta);

  std::vector<PetscReal> ni(mVtxCrd.size());

  VectorXd nodalPoints_x(mNumIntPnt);
  VectorXd nodalPoints_z(mNumIntPnt);

  int idx = 0;
  for (auto i = 0; i < mNumIntPtsEta; i++) {
    for (auto j = 0; j < mNumIntPtsEps; j++) {

      double eps = mIntCrdEps(j);
      double eta = mIntCrdEta(i);

      // reference mapping below uses [0,1]x[0,1] reference square
      double r = (eps + 1) / 2;
      double s = (eta + 1) / 2;

      // assumes right-hand rule vertex layout (same as PETSc)
      double v1x = mVtxCrd.row(0)[0];
      double v2x = mVtxCrd.row(0)[1];
      double v3x = mVtxCrd.row(0)[2];
      double v4x = mVtxCrd.row(0)[3];
      double v1z = mVtxCrd.row(1)[0];
      double v2z = mVtxCrd.row(1)[1];
      double v3z = mVtxCrd.row(1)[2];
      double v4z = mVtxCrd.row(1)[3];

      nodalPoints_x(idx) = v1x + (v2x - v1x) * r + (v4x + (v3x - v4x) * r - v1x - (v2x - v1x) * r) * s;
      nodalPoints_z(idx) = v1z + (v2z - v1z) * r + (v4z + (v3z - v4z) * r - v1z - (v2z - v1z) * r) * s;

      idx++;
    }
  }

  return std::make_tuple(nodalPoints_x, nodalPoints_z);
}


VectorXd Quad::interpolateLagrangePolynomials(
    const double eps, const double eta, const int p_order) {

  assert(p_order > 0 && p_order < 11);

  int n_points = (p_order + 1) * (p_order + 1);
  VectorXd gll_coeffs(n_points);
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
  Matrix2d inverse_Jacobian;
  double detJ;
  for (int i = 0; i < mNumIntPtsEta; i++) {
    for (int j = 0; j < mNumIntPtsEps; j++) {

      double eps = mIntCrdEps(j);
      double eta = mIntCrdEta(i);
      std::tie(inverse_Jacobian, detJ) = inverseJacobianAtPoint(eps, eta);
      val += field(j + i * mNumIntPtsEps) * mIntWgtEps(j) *
          mIntWgtEta(i) * detJ;

    }
  }

  return val;

}
void Quad::attachReceiver(std::vector<std::unique_ptr<Receiver>> &receivers) {

  for (auto &rec : receivers) { }


}


