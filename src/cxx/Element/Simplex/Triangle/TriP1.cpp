#include <vector>

#include <TriP1.h>

using namespace Eigen;

bool TriP1::checkHull(double x, double z,
                         const Ref<const Matrix<double,mNumVtx,mNumDim>> &vtx) {

  double x1 = vtx(0,0);
  double x2 = vtx(1,0);
  double x3 = vtx(2,0);
  double z1 = vtx(0,1);
  double z2 = vtx(1,1);
  double z3 = vtx(2,1);

  // see triangle_element_metrics.py
  // check barycentric coordinates of the point relative to this
  // triangle. If we violate barycentric coorindate assumptions
  // (l123 >= 0, l1+l2+l3 = 1) we are not inside the hull.
  double l1 = ((x - x3)*(z2 - z3) + (x2 - x3)*(-z + z3))/((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3));
  double l2 = (-x*(z1 - z3) + x3*(z1 - z3) + z*(x1 - x3) - z3*(x1 - x3)) /
    ((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3));
  double l3 = 1-l1-l2;
    
  if(l1 < 0) return false;
  if(l2 < 0) return false;
  if(l3 < 0) return false;
    
  if(fabs(1-(l1+l2+l3)) < 1e-8) return false;

  // no assumptions violated, inside triangle
  return true;
    
}

Vector2d TriP1::inverseCoordinateTransform(double x, double z,
                                              const Ref<const Matrix<double,mNumVtx,mNumDim>> &vtx) {

  double x1 = vtx(0,0);
  double x2 = vtx(1,0);
  double x3 = vtx(2,0);
  double z1 = vtx(0,1);
  double z2 = vtx(1,1);
  double z3 = vtx(2,1);
  
  
  // see triangle_element_metrics.py (and https://en.wikipedia.org/wiki/Barycentric_coordinate_system)
    
  double r  = ((-x + x3)*(z2 - z3)*((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)) + (x2 - x3)*(z - z3)*((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)) + (x*(z1 - z3) - x3*(z1 - z3) - z*(x1 - x3) + z3*(x1 - x3))*(x*(z1 - z3) - x3*(z1 - z3) - z*(x1 - x3) + z*(x2 - x3) + z3*(x1 - x3) - z3*(x2 - x3) + (-x + x3)*(z2 - z3) + (x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)))/pow((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3),2.0);
  double s = ((-x + x3)*(z2 - z3)*((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)) + (x2 - x3)*(z - z3)*((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)) + (x*(z1 - z3) - x3*(z1 - z3) - z*(x1 - x3) + z3*(x1 - x3))*(x*(z1 - z3) - x3*(z1 - z3) - z*(x1 - x3) + z*(x2 - x3) + z3*(x1 - x3) - z3*(x2 - x3) + (-x + x3)*(z2 - z3) + (x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)))/pow((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3),2.0);

  Vector2d solution {r, s};
  return solution;
}

std::tuple<Matrix2d,PetscReal> TriP1::inverseJacobian(const Ref<const Matrix<double,mNumVtx,mNumDim>> &vtx) {    
    
  double v1x = vtx(0,0);
  double v2x = vtx(1,0);
  double v3x = vtx(2,0);
                 
  double v1z = vtx(0,1);
  double v2z = vtx(1,1);
  double v3z = vtx(2,1);

  // tranform matrix for reference triangle (-1,-1),(1,-1),(-1,1)
  // See Hesthaven & Warburton pg. 172
  auto v2mv1_x = (v2x-v1x)/2;
  auto v2mv1_z = (v2z-v1z)/2;
  auto dxdr = v2mv1_x;
  auto dzdr = v2mv1_z;
  auto v3mv1_x = (v3x-v1x)/2;
  auto v3mv1_z = (v3z-v1z)/2;
  auto dxds = v3mv1_x;
  auto dzds = v3mv1_z;
    
  Matrix<double,2,2> J;
  J << dxdr,dzdr
    ,dxds,dzds;

  auto J_inv = 1/J.determinant();
  auto detJ = 1/J_inv;
    
  auto rx = J_inv*dzds;
  auto rz = -J_inv*dxds;
  auto sx = -J_inv*dzdr;
  auto sz = J_inv*dxdr;
    
  Matrix<double,2,2> inverseJacobian;
  inverseJacobian << rx, sx,
    rz, sz;
  return std::make_tuple(inverseJacobian, detJ);
    
}

Vector3d TriP1::interpolateAtPoint(double r, double s) {
  // interp_vector=[-r/2 - s/2, r/2 + 1/2, s/2 + 1/2]
  // barycentric coordinates l1,l2,l3 (l1+l2+l3=1) give us a linear weighting between coordinates (r,s). By computing the barycentric weights, we can compute the interpolated velocity = [l1,l2,l3].dot([v1,v2,v3]) where vN is the vertex velocity.
  // T*[l1;l2] = xy - v3, where xy is the desired point, and v3 is the vertex 3
  // thus l12 = T**(-1)*xy - T**(-1)*v3
  // For reference triangle
  // Tref = sym.Matrix([[0, 2],
  //                    [-2,-2]])
  // See triangle_element_metrics.py: `interp_vector`
  Vector3d interpolator;
  interpolator <<
    -r/2.0 - s/2.0, r/2.0 + 1.0/2.0, s/2.0 + 1.0/2.0;
  return interpolator;
}

// global x-z points on all nodes
std::tuple<VectorXd,VectorXd> TriP1::buildNodalPoints(const Ref<const VectorXd> &intCrdR,
                                                         const Ref<const VectorXd> &intCrdS,
                                                         const Ref<const Matrix<double,mNumVtx,mNumDim>> &vtx) {
		
  std::vector<PetscReal> ni(intCrdR.size());
	
  VectorXd nodalPoints_x(intCrdR.size());
  VectorXd nodalPoints_z(intCrdR.size());

  double x1 = vtx(0,0);
  double x2 = vtx(1,0);
  double x3 = vtx(2,0);
  double z1 = vtx(0,1);
  double z2 = vtx(1,1);
  double z3 = vtx(2,1);

  for(auto n = 0; n < intCrdR.size(); n++) {
    auto r = intCrdR[n];
    auto s = intCrdS[n];
    // map from reference triangle point (r,s) to this triangle (x,z)
    auto xn=-x1*(r + s)/2 + x2*(r + 1)/2 + x3*(s + 1)/2;
    auto zn=-z1*(r + s)/2 + z2*(r + 1)/2 + z3*(s + 1)/2;
    nodalPoints_x(n) = xn;
    nodalPoints_z(n) = zn;
  }    

  return std::make_tuple(nodalPoints_x,nodalPoints_z);
}
