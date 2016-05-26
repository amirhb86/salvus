#include <vector>
#include <HexP1.h>

using namespace Eigen;

bool HexP1::checkHull(const double x, const double y, const double z,
                      const Ref<const Matrix<double,mNumVtx,mNumDim>> &vtx)

{
  // assumes right-hand rule vertex layout (same as PETSc)
  double v0x = vtx(0,0);
  double v1x = vtx(1,0);
  double v2x = vtx(2,0);
  double v3x = vtx(3,0);
  double v4x = vtx(4,0);
  double v5x = vtx(5,0);
  double v6x = vtx(6,0);
  double v7x = vtx(7,0);
               
  double v0y = vtx(0,1);
  double v1y = vtx(1,1);
  double v2y = vtx(2,1);
  double v3y = vtx(3,1);
  double v4y = vtx(4,1);
  double v5y = vtx(5,1);
  double v6y = vtx(6,1);
  double v7y = vtx(7,1);
               
  double v0z = vtx(0,2);
  double v1z = vtx(1,2);
  double v2z = vtx(2,2);
  double v3z = vtx(3,2);
  double v4z = vtx(4,2);
  double v5z = vtx(5,2);
  double v6z = vtx(6,2);
  double v7z = vtx(7,2);

  // sides: bottom, top, front, back, left, right (all right-hand rule with normal outwards)
  std::vector<std::string> side_names = {
    "bottom",
    "top   ",
    "front ",
    "back  ",
    "left  ",
    "right "
  };
  std::vector<std::vector<double>> sides_x = {
    {v0x,v1x,v2x,v3x}, // bottom
    {v4x,v5x,v6x,v7x}, // top
    {v0x,v3x,v5x,v4x}, // front
    {v2x,v1x,v7x,v6x}, // back
    {v0x,v4x,v7x,v1x}, // left
    {v2x,v6x,v5x,v3x}  // right 
  };
  std::vector<std::vector<double>> sides_y = {
    {v0y,v1y,v2y,v3y}, // bottom
    {v4y,v5y,v6y,v7y}, // top
    {v0y,v3y,v5y,v4y}, // front
    {v2y,v1y,v7y,v6y}, // back
    {v0y,v4y,v7y,v1y}, // left
    {v2y,v6y,v5y,v3y} // right
  };

  std::vector<std::vector<double>> sides_z = {
    {v0z,v1z,v2z,v3z}, // bottom
    {v4z,v5z,v6z,v7z}, // top
    {v0z,v3z,v5z,v4z}, // front
    {v2z,v1z,v7z,v6z}, // back
    {v0z,v4z,v7z,v1z}, // left
    {v2z,v6z,v5z,v3z} // right
  };

  Vector3d test_point;
  test_point << x, y, z;
  
  for(int sideid=0;sideid<sides_z.size();sideid++) {
  
    Vector3d pt_a,pt_b,pt_d;
    pt_a << sides_x[sideid][0], sides_y[sideid][0], sides_z[sideid][0];
    pt_b << sides_x[sideid][1], sides_y[sideid][1], sides_z[sideid][1];
    pt_d << sides_x[sideid][3], sides_y[sideid][3], sides_z[sideid][3];
    auto line_1 = pt_b - pt_a;
    auto line_2 = pt_d - pt_a;
    auto normal_side = line_1.cross(line_2);
    auto normal_side_dot_test_point = normal_side.dot(test_point-pt_a);
    // if test point is colinear with the side's outward normal, it
    // will have a positive dot-product. This means that the point is
    // outside the hull of the hexahedra. This assumes planar sides.
    if(normal_side_dot_test_point > 0) {      
      return false;
    }
  }
  
  return true;
}

/** Derivation for reference to "real" element mapping
 *  reference hex, with r,s,t=[-1,1]x[-1,1]x[-1,1]
 *  
 *                                   (v7)                    (v6)
 *                                     /--------------d------/+
 *                                  /--|                  /-- |
 *                              /---   |              /---    |
 *                           /--       |       f   /--        |
 *                        /--          |         --           |
 *                 (v4) +---------------b--------+ (v5)       |
 *                      |              |         |            |
 *                      |              |         |            |
 *                      |              |       g |            |
 *    ^                 |              |         |            |
 *    | (t)             |              |         |            |
 *    |                 |            --+-------- |----c-------+ (v2)
 *    |                 |        ---/ (v1)       |         /--
 *    |                 |     --/              e |     /---
 *    |         (s)     | ---/                   |  /--
 *    |         /-      +/--------------a--------+--
 *    |     /---       (v0)                    (v3)
 *    | /---
 *    +-----------------> (r)
 * 
 *  
 * // note (r-1)/2.0 to bring range from [-1,1] to [0,1]
 * a = v0 + (v3-v0)*(r+1.0)/2.0
 * b = v4 + (v5-v4)*(r+1.0)/2.0
 * c = v1 + (v2-v1)*(r+1.0)/2.0
 * d = v7 + (v6-v7)*(r+1.0)/2.0
 * 
 * e = a + (c-a)*(s+1.0)/2.0
 * f = b + (d-b)*(s+1.0)/2.0
 * 
 * g = e + (f-e)*(t+1.0)/2.0
 * using sympy to expand `g` as a function of only vertices v0..7 and r,s,t gives the formula below
 */
PetscReal referenceToElementMapping(PetscReal v0, PetscReal v1, PetscReal v2, PetscReal v3, PetscReal v4, PetscReal v5, PetscReal v6, PetscReal v7, PetscReal r, PetscReal s, PetscReal t) {
  // represents `g` from above
  // return v0 - 0.5*(r + 1.0)*(v0 - v3) - 0.5*(s + 1.0)*(v0 - v1 - 0.5*(r + 1.0)*(v0 - v3) + 0.5*(r + 1.0)*(v1 - v2)) + 0.5*(t + 1.0)*(-v0 + v4 + 0.5*(r + 1.0)*(v0 - v3) - 0.5*(r + 1.0)*(v4 - v5) + 0.5*(s + 1.0)*(v0 - v1 - 0.5*(r + 1.0)*(v0 - v3) + 0.5*(r + 1.0)*(v1 - v2)) + 0.5*(s + 1.0)*(-v4 + v7 + 0.5*(r + 1.0)*(v4 - v5) + 0.5*(r + 1.0)*(v6 - v7)));
  return v0 + 0.5*(r + 1.0)*(-v0 + v3) + 0.5*(s + 1.0)*(-v0 + v1 - 0.5*(r + 1.0)*(-v0 + v3) + 0.5*(r + 1.0)*(-v1 + v2)) + 0.5*(t + 1.0)*(-v0 + v4 - 0.5*(r + 1.0)*(-v0 + v3) + 0.5*(r + 1.0)*(-v4 + v5) - 0.5*(s + 1.0)*(-v0 + v1 - 0.5*(r + 1.0)*(-v0 + v3) + 0.5*(r + 1.0)*(-v1 + v2)) + 0.5*(s + 1.0)*(-v4 + v7 - 0.5*(r + 1.0)*(-v4 + v5) + 0.5*(r + 1.0)*(v6 - v7)));
}

Vector3d HexP1::inverseCoordinateTransform(const double x_real, const double y_real, const double z_real,
                                           const Ref<const Matrix<double,mNumVtx,mNumDim>> &vtx) {

  // assumes right-hand rule vertex layout (same as PETSc)
  double v0x = vtx(0,0);
  double v1x = vtx(1,0);
  double v2x = vtx(2,0);
  double v3x = vtx(3,0);
  double v4x = vtx(4,0);
  double v5x = vtx(5,0);
  double v6x = vtx(6,0);
  double v7x = vtx(7,0);
               
  double v0y = vtx(0,1);
  double v1y = vtx(1,1);
  double v2y = vtx(2,1);
  double v3y = vtx(3,1);
  double v4y = vtx(4,1);
  double v5y = vtx(5,1);
  double v6y = vtx(6,1);
  double v7y = vtx(7,1);
               
  double v0z = vtx(0,2);
  double v1z = vtx(1,2);
  double v2z = vtx(2,2);
  double v3z = vtx(3,2);
  double v4z = vtx(4,2);
  double v5z = vtx(5,2);
  double v6z = vtx(6,2);
  double v7z = vtx(7,2);

  // Using Newton iterations
  // https://en.wikipedia.org/wiki/Newton%27s_method#Nonlinear_systems_of_equations
  // J_F(xn)(x_{n+1} - x_n) = -F(x_n)
  // Solve for x_{n+1}
  // where J_F(x_n) is jacobian:
  // https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
  double tol = 1e-6;
  int num_iter = 0;
  Vector3d solution{0.0, 0.0, 0.0}; // initial guess at (0.0,0.0)
  while (true) {

    double r = solution(0);
    double s = solution(1);
    double t = solution(2);

    Matrix3d jacobian_inverse_t;
    PetscReal detJ;

    // mapping from reference hex [-1,1]x[-1,1]x[-1,1] to *this* element
    PetscReal Tx = referenceToElementMapping(v0x,v1x,v2x,v3x,v4x,v5x,v6x,v7x,r,s,t);
    PetscReal Ty = referenceToElementMapping(v0y,v1y,v2y,v3y,v4y,v5y,v6y,v7y,r,s,t);
    PetscReal Tz = referenceToElementMapping(v0z,v1z,v2z,v3z,v4z,v5z,v6z,v7z,r,s,t);
     
    Vector3d objective_function{x_real - Tx, y_real - Ty, z_real - Tz};
    
    std::tie(jacobian_inverse_t,detJ) = inverseJacobianAtPoint(r,s,t,vtx);
    
    if ((objective_function.array().abs() < tol).all()) {
      return solution;
    } else {
      solution += (jacobian_inverse_t.transpose() * objective_function);
    }
    if (num_iter > 10) {
      std::cerr << "inverseCoordinateTransform: TOO MANY ITERATIONS!\n";
      exit(1);
    }
    num_iter++;

  }

}

std::tuple<Matrix3d, PetscReal> HexP1::inverseJacobianAtPoint(PetscReal r,
                                                              PetscReal s,
                                                              PetscReal t,
                                                              const Ref<const Matrix<double,mNumVtx,mNumDim>> &vtx) {

  // assumes right-hand rule vertex layout (same as PETSc)
  double v0x = vtx(0,0);
  double v1x = vtx(1,0);
  double v2x = vtx(2,0);
  double v3x = vtx(3,0);
  double v4x = vtx(4,0);
  double v5x = vtx(5,0);
  double v6x = vtx(6,0);
  double v7x = vtx(7,0);

  double v0y = vtx(0,1);
  double v1y = vtx(1,1);
  double v2y = vtx(2,1);
  double v3y = vtx(3,1);
  double v4y = vtx(4,1);
  double v5y = vtx(5,1);
  double v6y = vtx(6,1);
  double v7y = vtx(7,1);

  double v0z = vtx(0,2);
  double v1z = vtx(1,2);
  double v2z = vtx(2,2);
  double v3z = vtx(3,2);
  double v4z = vtx(4,2);
  double v5z = vtx(5,2);
  double v6z = vtx(6,2);
  double v7z = vtx(7,2);

  PetscReal dxdr=-0.5*v0x + 0.5*v3x + 0.5*(s + 1.0)*(0.5*v0x - 0.5*v1x + 0.5*v2x - 0.5*v3x) + 0.5*(t + 1.0)*(0.5*v0x - 0.5*v3x - 0.5*v4x + 0.5*v5x - 0.5*(s + 1.0)*(0.5*v0x - 0.5*v1x + 0.5*v2x - 0.5*v3x) + 0.5*(s + 1.0)*(0.5*v4x - 0.5*v5x + 0.5*v6x - 0.5*v7x));
PetscReal dydr=-0.5*v0y + 0.5*v3y + 0.5*(s + 1.0)*(0.5*v0y - 0.5*v1y + 0.5*v2y - 0.5*v3y) + 0.5*(t + 1.0)*(0.5*v0y - 0.5*v3y - 0.5*v4y + 0.5*v5y - 0.5*(s + 1.0)*(0.5*v0y - 0.5*v1y + 0.5*v2y - 0.5*v3y) + 0.5*(s + 1.0)*(0.5*v4y - 0.5*v5y + 0.5*v6y - 0.5*v7y));
PetscReal dzdr=-0.5*v0z + 0.5*v3z + 0.5*(s + 1.0)*(0.5*v0z - 0.5*v1z + 0.5*v2z - 0.5*v3z) + 0.5*(t + 1.0)*(0.5*v0z - 0.5*v3z - 0.5*v4z + 0.5*v5z - 0.5*(s + 1.0)*(0.5*v0z - 0.5*v1z + 0.5*v2z - 0.5*v3z) + 0.5*(s + 1.0)*(0.5*v4z - 0.5*v5z + 0.5*v6z - 0.5*v7z));
PetscReal dxds=-0.5*v0x + 0.5*v1x - 0.25*(r + 1.0)*(-v0x + v3x) + 0.25*(r + 1.0)*(-v1x + v2x) + 0.5*(t + 1.0)*(0.5*v0x - 0.5*v1x - 0.5*v4x + 0.5*v7x + 0.25*(r + 1.0)*(-v0x + v3x) - 0.25*(r + 1.0)*(-v1x + v2x) - 0.25*(r + 1.0)*(-v4x + v5x) + 0.25*(r + 1.0)*(v6x - v7x));
PetscReal dyds=-0.5*v0y + 0.5*v1y - 0.25*(r + 1.0)*(-v0y + v3y) + 0.25*(r + 1.0)*(-v1y + v2y) + 0.5*(t + 1.0)*(0.5*v0y - 0.5*v1y - 0.5*v4y + 0.5*v7y + 0.25*(r + 1.0)*(-v0y + v3y) - 0.25*(r + 1.0)*(-v1y + v2y) - 0.25*(r + 1.0)*(-v4y + v5y) + 0.25*(r + 1.0)*(v6y - v7y));
PetscReal dzds=-0.5*v0z + 0.5*v1z - 0.25*(r + 1.0)*(-v0z + v3z) + 0.25*(r + 1.0)*(-v1z + v2z) + 0.5*(t + 1.0)*(0.5*v0z - 0.5*v1z - 0.5*v4z + 0.5*v7z + 0.25*(r + 1.0)*(-v0z + v3z) - 0.25*(r + 1.0)*(-v1z + v2z) - 0.25*(r + 1.0)*(-v4z + v5z) + 0.25*(r + 1.0)*(v6z - v7z));
PetscReal dxdt=-0.5*v0x + 0.5*v4x - 0.25*(r + 1.0)*(-v0x + v3x) + 0.25*(r + 1.0)*(-v4x + v5x) - 0.25*(s + 1.0)*(-v0x + v1x - 0.5*(r + 1.0)*(-v0x + v3x) + 0.5*(r + 1.0)*(-v1x + v2x)) + 0.25*(s + 1.0)*(-v4x + v7x - 0.5*(r + 1.0)*(-v4x + v5x) + 0.5*(r + 1.0)*(v6x - v7x));
PetscReal dydt=-0.5*v0y + 0.5*v4y - 0.25*(r + 1.0)*(-v0y + v3y) + 0.25*(r + 1.0)*(-v4y + v5y) - 0.25*(s + 1.0)*(-v0y + v1y - 0.5*(r + 1.0)*(-v0y + v3y) + 0.5*(r + 1.0)*(-v1y + v2y)) + 0.25*(s + 1.0)*(-v4y + v7y - 0.5*(r + 1.0)*(-v4y + v5y) + 0.5*(r + 1.0)*(v6y - v7y));
PetscReal dzdt=-0.5*v0z + 0.5*v4z - 0.25*(r + 1.0)*(-v0z + v3z) + 0.25*(r + 1.0)*(-v4z + v5z) - 0.25*(s + 1.0)*(-v0z + v1z - 0.5*(r + 1.0)*(-v0z + v3z) + 0.5*(r + 1.0)*(-v1z + v2z)) + 0.25*(s + 1.0)*(-v4z + v7z - 0.5*(r + 1.0)*(-v4z + v5z) + 0.5*(r + 1.0)*(v6z - v7z));
    
  Matrix3d J;
  J <<
    dxdr, dydr, dzdr,
    dxds, dyds, dzds,
    dxdt, dydt, dzdt;

  PetscReal detJ = J.determinant();

  // detJ * invJ_ = [ dyds*dzdt - dydt*dzds, -dydr*dzdt + dydt*dzdr,  dydr*dzds - dyds*dzdr]
  //                [-dxds*dzdt + dxdt*dzds,  dxdr*dzdt - dxdt*dzdr, -dxdr*dzds + dxds*dzdr]
  //                [ dxds*dydt - dxdt*dyds, -dxdr*dydt + dxdt*dydr,  dxdr*dyds - dxds*dydr]  
  
  // need to test if this quasi-analytical formula is faster than simply J.inverse();
  
  double rx = (1 / detJ) * (dyds*dzdt - dydt*dzds);
  double sx = (1 / detJ) * (-dydr*dzdt + dydt*dzdr);
  double tx = (1 / detJ) * (dydr*dzds - dyds*dzdr);

  double ry = (1 / detJ) * (-dxds*dzdt + dxdt*dzds);
  double sy = (1 / detJ) * ( dxdr*dzdt - dxdt*dzdr);
  double ty = (1 / detJ) * (-dxdr*dzds + dxds*dzdr);
  
  double rz = (1 / detJ) * (dxds*dydt - dxdt*dyds);
  double sz = (1 / detJ) * (-dxdr*dydt + dxdt*dydr);
  double tz = (1 / detJ) * (dxdr*dyds - dxds*dydr);

  Matrix3d inverseJacobian;
  inverseJacobian <<
    rx, sx, tx,
    ry, sy, ty,
    rz, sz, tz;

  return std::make_tuple(inverseJacobian, detJ);
}

Eigen::VectorXd HexP1::interpolateAtPoint(double r, double s, double t) {
  /**  
   *  reference hex, with r,s,t=[-1,1]x[-1,1]x[-1,1]
   * 
   *                                      (v7)                    (v6)
   *                                        /--------------d------/+
   *                                     /--|                  /-- |
   *                                 /---   |              /---    |
   *                              /--       |       f   /--        |
   *                           /--          |         --           |
   *                    (v4) +---------------b--------+ (v5)       |        \gamma
   *                         |              |         |            |          ^
   *                         |              |         |            |          |
   *                         |              |       g |            |          |
   *    ^                    |              |         |            |          |
   *    | (t)                |              |         |            |          |
   *    |                    |            --+-------- |----c-------+ (v2)     |
   *    |                    |        ---/ (v1)       |         /--          ---
   *    |                    |     --/              e |     /---               /> \beta
   *    |         (s)        | ---/                   |  /--                 /-
   *    |         /-         +/--------------a--------+--                  --
   *    |     /---          (v0)                    (v3)
   *    | /---               |--------------> \alpha
   *    +-----------------> (r)
   * 
   * alpha = (r+1)/2.0
   * beta = (s+1)/2.0
   * gamma = (t+1)/2.0
   * 
   * a = alpha*v3 + (1-alpha)*v0
   * b = alpha*v5 + (1-alpha)*v4
   * c = alpha*v2 + (1-alpha)*v1
   * d = alpha*v6 + (1-alpha)*v7
   * 
   * e = beta*c + (1-beta)*a
   * f = beta*d + (1-beta)*b
   * 
   * g = gamma*f + (1-gamma)*e
   
   * group vf by terms v0..v7                                       
   * = [v0,...,v7].dot([-0.125*r*s*t + 0.125*r*s + ... - 0.125*s - 0.125*t + 0.125,
   *                    ...,
   *                    -0.125*r*s*t - 0.125*r*s - ... + 0.125*s + 0.125*t + 0.125]
   */
  
  Eigen::VectorXd interpolator(8);
  interpolator <<
    -0.125*r*s*t + 0.125*r*s + 0.125*r*t - 0.125*r + 0.125*s*t - 0.125*s - 0.125*t + 0.125,
    0.125*r*s*t - 0.125*r*s + 0.125*r*t - 0.125*r - 0.125*s*t + 0.125*s - 0.125*t + 0.125,
    -0.125*r*s*t + 0.125*r*s - 0.125*r*t + 0.125*r - 0.125*s*t + 0.125*s - 0.125*t + 0.125,
    0.125*r*s*t - 0.125*r*s - 0.125*r*t + 0.125*r + 0.125*s*t - 0.125*s - 0.125*t + 0.125,
    0.125*r*s*t + 0.125*r*s - 0.125*r*t - 0.125*r - 0.125*s*t - 0.125*s + 0.125*t + 0.125,
    -0.125*r*s*t - 0.125*r*s + 0.125*r*t + 0.125*r - 0.125*s*t - 0.125*s + 0.125*t + 0.125,
    0.125*r*s*t + 0.125*r*s + 0.125*r*t + 0.125*r + 0.125*s*t + 0.125*s + 0.125*t + 0.125,
    -0.125*r*s*t - 0.125*r*s - 0.125*r*t - 0.125*r + 0.125*s*t + 0.125*s + 0.125*t + 0.125;
  return interpolator;
}

std::tuple<Eigen::VectorXd,
           Eigen::VectorXd,
           Eigen::VectorXd> HexP1::buildNodalPoints(const Eigen::Ref<const Eigen::VectorXd> &intCrdR,
                                                    const Eigen::Ref<const Eigen::VectorXd> &intCrdS,
                                                    const Eigen::Ref<const Eigen::VectorXd> &intCrdT,
                                                    const Eigen::Ref<const Eigen::Matrix<double,
                                                    mNumVtx,
                                                    mNumDim>> &vtx) {
  
  int num_pts_r = intCrdR.size();
  int num_pts_s = intCrdS.size();
  int num_pts_t = intCrdT.size();
  int num_pts = num_pts_r * num_pts_s * num_pts_t;
  VectorXd nodalPoints_x(num_pts);
  VectorXd nodalPoints_y(num_pts);
  VectorXd nodalPoints_z(num_pts);

  // assumes right-hand rule vertex layout (same as PETSc)
  double v0x = vtx(0,0);
  double v1x = vtx(1,0);
  double v2x = vtx(2,0);
  double v3x = vtx(3,0);
  double v4x = vtx(4,0);
  double v5x = vtx(5,0);
  double v6x = vtx(6,0);
  double v7x = vtx(7,0);

  double v0y = vtx(0,1);
  double v1y = vtx(1,1);
  double v2y = vtx(2,1);
  double v3y = vtx(3,1);
  double v4y = vtx(4,1);
  double v5y = vtx(5,1);
  double v6y = vtx(6,1);
  double v7y = vtx(7,1);

  double v0z = vtx(0,2);
  double v1z = vtx(1,2);
  double v2z = vtx(2,2);
  double v3z = vtx(3,2);
  double v4z = vtx(4,2);
  double v5z = vtx(5,2);
  double v6z = vtx(6,2);
  double v7z = vtx(7,2);
  
  int idx = 0;
  for (auto k = 0; k < num_pts_t; k++) {
    for (auto j = 0; j < num_pts_t; j++) {
      for (auto i = 0; i < num_pts_r; i++) {
        double r = intCrdR(i);
        double s = intCrdS(j);
        double t = intCrdT(k);

        // this represents `g` above (for each x,y,z) (copy paste from hexahedra_element_metrics.py)
        nodalPoints_x(idx) = referenceToElementMapping(v0x,v1x,v2x,v3x,v4x,v5x,v6x,v7x,r,s,t);
        nodalPoints_y(idx) = referenceToElementMapping(v0y,v1y,v2y,v3y,v4y,v5y,v6y,v7y,r,s,t);
        nodalPoints_z(idx) = referenceToElementMapping(v0z,v1z,v2z,v3z,v4z,v5z,v6z,v7z,r,s,t);

        idx++;
      }
    }
  }
  return std::make_tuple(nodalPoints_x,nodalPoints_y,nodalPoints_z);
}
