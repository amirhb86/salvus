#include <vector>
#include <Eigen/Core>
#include <Element/HyperCube/QuadP1.h>
#include <stdexcept>

using namespace Eigen;

bool QuadP1::checkHull(const PetscReal x, const PetscReal y, const Ref<const QuadVtx>& vtx) {

  EIGEN_ASM_COMMENT("HI");
  /* As a bonus, check to see if the right hand rule is followed. */
  PetscReal detJac;
  RealMat2x2 invJac;
  inverseJacobianAtPoint(0, 0, vtx, detJac, invJac);
  if (detJac <= 0) {
    throw std::runtime_error("Poorly behaved TensorQuad detected by QuadP1");
  }

  PetscInt n_neg = 0;
  PetscInt n_pos = 0;
  std::vector<int> edge_mapping{0, 1, 2, 3, 0};
  RealVec2 test_point;
  test_point << x, y;
  for (PetscInt i = 0; i < mNumVtx; i++) {
    RealVec2 p0 = vtx.row(edge_mapping[i + 0]);
    RealVec2 p1 = vtx.row(edge_mapping[i + 1]);
    EIGEN_ASM_COMMENT("VECTORIZE BEGIN");
    RealVec2 v_seg = p1 - p0;
    EIGEN_ASM_COMMENT("VECTORIZE END");
    RealVec2 p_seg = test_point - p0;
    PetscReal x_0 = v_seg(0) * p_seg(1) - v_seg(1) * p_seg(0);
    if (x_0 < 0) {
      n_neg++;
    } else if (x_0 > 0) {
      n_pos++;
    } else { // if x_0 == 0.
      n_pos++; n_neg++;
    }
  }
  return n_neg == mNumVtx || n_pos == mNumVtx;
}

RealVec2 QuadP1::inverseCoordinateTransform(const PetscReal x,
                                            const PetscReal y,
                                            const Ref<const QuadVtx> &vtx) {

  PetscReal v1x = vtx(0,0);
  PetscReal v2x = vtx(1,0);
  PetscReal v3x = vtx(2,0);
  PetscReal v4x = vtx(3,0);
  PetscReal v1z = vtx(0,1);
  PetscReal v2z = vtx(1,1);
  PetscReal v3z = vtx(2,1);
  PetscReal v4z = vtx(3,1);

  // Using Newton iterations
  // https://en.wikipedia.org/wiki/Newton%27s_method#Nonlinear_systems_of_equations
  PetscReal tol = 1e-6;
  PetscInt num_iter = 0;
  RealVec2 solution {0.0, 0.0}; // initial guess at (0.0,0.0)

  while (true) {

    PetscReal r = solution(0);
    PetscReal s = solution(1);

    RealMat2x2 jacobian;

    // mapping from reference quad [-1,1]x[-1,1] to *this* element
    PetscReal Tx = v1x + (v2x - v1x) * (r + 1) / 2
        + (v4x + (v3x - v4x) * (r + 1) / 2 - v1x - (v2x - v1x) * (r + 1) / 2) * (s + 1) / 2;
    PetscReal Tz = v1z + (v2z - v1z) * (r + 1) / 2
        + (v4z + (v3z - v4z) * (r + 1) / 2 - v1z - (v2z - v1z) * (r + 1) / 2) * (s + 1) / 2;
    RealVec2 objective_function{x - Tx, y - Tz};

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
      throw std::runtime_error("inverseCoordinateTransform in QuadP1 failed to converge after "
                                   "10 iteration.");
    }
    num_iter++;

  }

}

//void QuadP1::inverseJacobianAtPoint(
//    PetscReal r, PetscReal s, const Ref<const QuadVtx> &vtx,
//    PetscReal &detJac, Eigen::Ref<RealMat2x2> invJac) {
//
//  /*
//   * Below I've been experimenting with some optimizations by looking at the assembly. It
//   * looks like using Eigen, along with these static functions, cuts the number of
//   * assembly instructions by half. Still, I can't seem to get it to vectorize properly.
//   * May be a hardware limitation?
//   * Remeber: to compile to assembly: g++ -std=c++11 -c -S -O2.
//   */
//  EIGEN_ASM_COMMENT("JACOBIAN NEW");
//  RealMat2x2 jac;
//  Matrix<PetscReal,2,4> mult;
//  jac = (mult << dn0dr(r), dn1dr(r), dn2dr(r), dn3dr(r),
//                 dn0ds(s), dn1ds(s), dn2ds(s), dn3ds(s)).finished() * vtx;
//  detJac = jac.determinant();
//  invJac = jac.inverse();
//  EIGEN_ASM_COMMENT("END_JACOBIAN_NEW");
//
////  EIGEN_ASM_COMMENT("JACOBIAN OLD");
////  PetscReal v1x = vtx(0, 0);
////  PetscReal v2x = vtx(1, 0);
////  PetscReal v3x = vtx(2, 0);
////  PetscReal v4x = vtx(3, 0);
////  PetscReal v1y = vtx(0, 1);
////  PetscReal v2y = vtx(1, 1);
////  PetscReal v3y = vtx(2, 1);
////  PetscReal v4y = vtx(3, 1);
////
////  PetscReal r_ = (r + 1.0);
////  PetscReal s_ = (s + 1.0);
////  detJac = -r_ * v1x * v3y / 8 + r_ * v1x * v4y / 8 + r_ * v1y * v3x / 8 - r_ * v1y *
////      v4x / 8 + r_ * v2x * v3y / 8 - r_ * v2x * v4y / 8 - r_ * v2y * v3x / 8 + r_ * v2y * v4x /
////      8 - s_ * v1x * v2y / 8 + s_ * v1x * v3y / 8 + s_ * v1y * v2x / 8 - s_ * v1y * v3x / 8 - s_ *
////      v2x * v4y / 8 + s_ * v2y * v4x / 8 + s_ * v3x * v4y / 8 - s_ * v3y * v4x / 8 + v1x * v2y / 4
////      - v1x * v4y / 4 - v1y * v2x / 4 + v1y * v4x / 4 + v2x * v4y / 4 - v2y * v4x / 4;
////
////  // J      = [[dx/dr, dy/dr], [dx/ds, dy/ds]]
////  // J^{-1} = [[dr/dx, ds/dx], [dr/dz, ds/dz]]
////  PetscReal rx = (1 / detJac) * (r_ * (v1y - v2y) / 4 + r_ * (v3y - v4y) / 4 - v1y / 2 + v4y / 2);
////  PetscReal rz = (1 / detJac) * (-r_ * (v1x - v2x) / 4 - r_ * (v3x - v4x) / 4 + v1x / 2 - v4x / 2);
////  PetscReal sx = (1 / detJac) * (-s_ * (v1y - v2y + v3y - v4y) / 4 + v1y / 2 - v2y / 2);
////  PetscReal sz = (1 / detJac) * (s_ * (v1x - v2x + v3x - v4x) / 4 - v1x / 2 + v2x / 2);
////
////  invJac <<  rx, sx, rz, sz;
////  std::cout <<invJac << std::endl;
////  EIGEN_ASM_COMMENT("END JACOBIAN OLD");
//
//}

//RealVec4 QuadP1::interpolateAtPoint(const PetscReal r, const PetscReal s) {
//
//  // see element_metrics.py
//  //  interpolation for quadrilaterals based on right-hand rule quad with
//  // velocity v0,v1,v2,v3 and r(x) s(y,or,z)
//  //
//  //  v3               vb          v2
//  //  +-----------------+----------+ \s
//  //  |                 |          |   ^       \alpha = (r-(-1))/2
//  //  |                 |          |   |       \bs  = (s-(-1))/2
//  //  |              vf + \bs    |   |         \va = alpha(v1) + (1-alpha)v0
//  //  |                 |  ^       |   |       \vb = alpha(v2) + (1-alpha)v3
//  //  |                 |  |       |           \vf = bs va + (1-bs) vb
//  //  |                 |  |       |           group vf by terms v0..v3
//  //  |                 |  |       |           = [v0,v1,v2,v3].dot([-0.25*r*s - 0.25*r + 0.25*s + 0.25,
//  //  |              va |  |       |                                0.25*r*s + 0.25*r + 0.25*s + 0.25,
//  //  +-----------------+----------+                                -0.25*r*s + 0.25*r - 0.25*s + 0.25,
//  //  v0 -------------->\alpha    v1  --->\r                      0.25*r*s - 0.25*r - 0.25*s + 0.25]
//  // ----------------------------------------------------------------------------------------------------------
//
//  RealVec4 interpolator;
//  interpolator <<
//      +0.25*r*s - 0.25*r - 0.25*s + 0.25,
//      -0.25*r*s + 0.25*r - 0.25*s + 0.25,
//      +0.25*r*s + 0.25*r + 0.25*s + 0.25,
//      -0.25*r*s - 0.25*r + 0.25*s + 0.25;
//  return interpolator;
//
//}


std::tuple<RealVec,RealVec> QuadP1::buildNodalPoints(const Ref<const RealVec> &intCrdR,
                                                     const Ref<const RealVec> &intCrdS,
                                                     const Ref<const QuadVtx>& vtx) {

  PetscInt num_pnt_r = intCrdR.size();
  PetscInt num_pnt_s = intCrdS.size();
  PetscInt num_pnt = num_pnt_r * num_pnt_s;
  RealVec nodalPoints_x(num_pnt);
  RealVec nodalPoints_z(num_pnt);

  // assumes right-hand rule vertex layout (same as PETSc)
  PetscReal v1x = vtx(0, 0);
  PetscReal v2x = vtx(1, 0);
  PetscReal v3x = vtx(2, 0);
  PetscReal v4x = vtx(3, 0);
  PetscReal v1z = vtx(0, 1);
  PetscReal v2z = vtx(1, 1);
  PetscReal v3z = vtx(2, 1);
  PetscReal v4z = vtx(3, 1);

  PetscInt idx = 0;
  for (auto i = 0; i < num_pnt_s; i++) {
    for (auto j = 0; j < num_pnt_r; j++) {

      PetscReal eps = intCrdR(j);
      PetscReal eta = intCrdS(i);

      // reference mapping below uses [0,1]x[0,1] reference square
      PetscReal r = (eps + 1) / 2;
      PetscReal s = (eta + 1) / 2;
      
      nodalPoints_x(idx) = v1x + (v2x - v1x) * r + (v4x + (v3x - v4x) * r - v1x - (v2x - v1x) * r) * s;
      nodalPoints_z(idx) = v1z + (v2z - v1z) * r + (v4z + (v3z - v4z) * r - v1z - (v2z - v1z) * r) * s;

      idx++;
    }
  }

  return std::make_tuple(nodalPoints_x, nodalPoints_z);

}


