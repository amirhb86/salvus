#include <vector>
#include <Element/HyperCube/QuadP1.h>

using namespace Eigen;

bool QuadP1::checkHull(const double x, const double y,
                       const Ref<const Matrix<double,mNumVtx,mNumDim>> &vtx) {
  int n_neg = 0;
  int n_pos = 0;
  std::vector<int> edge_mapping{0, 1, 2, 3, 0};
  Vector2d test_point;
  test_point << x, y;
  for (int i = 0; i < mNumVtx; i++) {
    Vector2d p0 = vtx.row(edge_mapping[i + 0]);
    Vector2d p1 = vtx.row(edge_mapping[i + 1]);
    Vector2d v_seg = p1 - p0;
    Vector2d p_seg = test_point - p0;
    double x_0 = v_seg(0) * p_seg(1) - v_seg(1) * p_seg(0);
    if (x_0 < 0) {
      n_neg++;
    } else if (x_0 > 0) {
      n_pos++;
    } else {
      n_pos++; n_neg++;
    }
  }
  return n_neg == mNumVtx || n_pos == mNumVtx;
}

Eigen::Vector2d QuadP1::inverseCoordinateTransform(const double x,
                                                   const double y,
                                                   const Ref<const Matrix<double,mNumVtx,mNumDim>> &vtx) {

  double v1x = vtx(0,0);
  double v2x = vtx(1,0);
  double v3x = vtx(2,0);
  double v4x = vtx(3,0);
  double v1z = vtx(0,1);
  double v2z = vtx(1,1);
  double v3z = vtx(2,1);
  double v4z = vtx(3,1);

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
    Vector2d objective_function{x - Tx, y - Tz};

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

std::tuple<Matrix2d, PetscReal> QuadP1::inverseJacobianAtPoint(
    PetscReal r, PetscReal s, const Ref<const Matrix<double,mNumVtx,mNumDim>> &vtx) {

  double v1x = vtx(0, 0);
  double v2x = vtx(1, 0);
  double v3x = vtx(2, 0);
  double v4x = vtx(3, 0);
  double v1z = vtx(0, 1);
  double v2z = vtx(1, 1);
  double v3z = vtx(2, 1);
  double v4z = vtx(3, 1);

  double r_ = (r + 1.0);
  double s_ = (s + 1.0);
  double detJ = -r_ * v1x * v3z / 8 + r_ * v1x * v4z / 8 + r_ * v1z * v3x / 8 - r_ * v1z * v4x / 8 + r_ * v2x * v3z / 8 -
      r_ * v2x * v4z / 8 - r_ * v2z * v3x / 8 + r_ * v2z * v4x / 8 - s_ * v1x * v2z / 8 + s_ * v1x * v3z / 8 +
      s_ * v1z * v2x / 8 - s_ * v1z * v3x / 8 - s_ * v2x * v4z / 8 + s_ * v2z * v4x / 8 + s_ * v3x * v4z / 8 -
      s_ * v3z * v4x / 8 + v1x * v2z / 4 - v1x * v4z / 4 - v1z * v2x / 4 + v1z * v4x / 4 + v2x * v4z / 4 -
      v2z * v4x / 4;

  // J        =   [dx/dr, dy/dr;
  //               dx/ds, dy/ds]
  // J^{-1}   =   [dr/dx ds/dx;
  //               dr/dz ds/dz]
  double rx = (1 / detJ) * (r_ * (v1z - v2z) / 4 + r_ * (v3z - v4z) / 4 - v1z / 2 + v4z / 2);
  double rz = (1 / detJ) * (-r_ * (v1x - v2x) / 4 - r_ * (v3x - v4x) / 4 + v1x / 2 - v4x / 2);
  double sx = (1 / detJ) * (-s_ * (v1z - v2z + v3z - v4z) / 4 + v1z / 2 - v2z / 2);
  double sz = (1 / detJ) * (s_ * (v1x - v2x + v3x - v4x) / 4 - v1x / 2 + v2x / 2);

  Matrix2d inverseJacobian;
  inverseJacobian <<  rx, sx,
      rz, sz;

  return std::make_tuple(inverseJacobian, detJ);
}

Eigen::Vector4d QuadP1::interpolateAtPoint(const double r, const double s) {

  // see element_metrics.py
  //  interpolation for quadrilaterals based on right-hand rule quad with
  // velocity v0,v1,v2,v3 and r(x) s(y,or,z)
  //
  //  v3               vb          v2
  //  +-----------------+----------+ \s
  //  |                 |          |   ^       \alpha = (r-(-1))/2
  //  |                 |          |   |       \bs  = (s-(-1))/2
  //  |              vf + \bs    |   |       va = alpha(v1) + (1-alpha)v0
  //  |                 |  ^       |   |       vb = alpha(v2) + (1-alpha)v3
  //  |                 |  |       |           vf = bs va + (1-bs) vb
  //  |                 |  |       |           group vf by terms v0..v3
  //  |                 |  |       |           = [v0,v1,v2,v3].dot([-0.25*r*s - 0.25*r + 0.25*s + 0.25,
  //  |              va |  |       |                                0.25*r*s + 0.25*r + 0.25*s + 0.25,
  //  +-----------------+----------+                                -0.25*r*s + 0.25*r - 0.25*s + 0.25,
  //  v0 -------------->\alpha    v1  --->\r                      0.25*r*s - 0.25*r - 0.25*s + 0.25]
  // ----------------------------------------------------------------------------------------------------------
  Eigen::Vector4d interpolator;
  interpolator <<
      0.25*r*s - 0.25*r - 0.25*s + 0.25,
      -0.25*r*s + 0.25*r - 0.25*s + 0.25,
      0.25*r*s + 0.25*r + 0.25*s + 0.25,
      -0.25*r*s - 0.25*r + 0.25*s + 0.25;
  return interpolator;

}


std::tuple<VectorXd,VectorXd> QuadP1::buildNodalPoints(const Eigen::Ref<const Eigen::VectorXd> &intCrdR,
                                                       const Eigen::Ref<const Eigen::VectorXd> &intCrdS,
                                                       const Eigen::Ref<const Eigen::Matrix<double,
                                                                                            mNumVtx,
                                                                                            mNumDim>> &vtx) {

  int num_pnt_r = intCrdR.size();
  int num_pnt_s = intCrdS.size();
  int num_pnt = num_pnt_r * num_pnt_s;
  VectorXd nodalPoints_x(num_pnt);
  VectorXd nodalPoints_z(num_pnt);

  // assumes right-hand rule vertex layout (same as PETSc)
  double v1x = vtx(0, 0);
  double v2x = vtx(1, 0);
  double v3x = vtx(2, 0);
  double v4x = vtx(3, 0);
  double v1z = vtx(0, 1);
  double v2z = vtx(1, 1);
  double v3z = vtx(2, 1);
  double v4z = vtx(3, 1);
  
  int idx = 0;
  for (auto i = 0; i < num_pnt_s; i++) {
    for (auto j = 0; j < num_pnt_r; j++) {

      double eps = intCrdR(j);
      double eta = intCrdS(i);

      // reference mapping below uses [0,1]x[0,1] reference square
      double r = (eps + 1) / 2;
      double s = (eta + 1) / 2;
      
      nodalPoints_x(idx) = v1x + (v2x - v1x) * r + (v4x + (v3x - v4x) * r - v1x - (v2x - v1x) * r) * s;
      nodalPoints_z(idx) = v1z + (v2z - v1z) * r + (v4z + (v3z - v4z) * r - v1z - (v2z - v1z) * r) * s;

      idx++;
    }
  }

  return std::make_tuple(nodalPoints_x, nodalPoints_z);

}


