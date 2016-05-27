#include <vector>
#include <Element/Simplex/TetP1.h>

using namespace Eigen;

bool TetP1::checkHull(double x, double y, double z,
                      const Ref<const Matrix<double,mNumVtx,mNumDim>> &vtx) {

  auto x1 = vtx(0,0);
  auto x2 = vtx(1,0);
  auto x3 = vtx(2,0);
  auto x4 = vtx(3,0);
            
  auto y1 = vtx(0,1);
  auto y2 = vtx(1,1);
  auto y3 = vtx(2,1);
  auto y4 = vtx(3,1);
            
  auto z1 = vtx(0,2);
  auto z2 = vtx(1,2);
  auto z3 = vtx(2,2);
  auto z4 = vtx(3,2);

  // see tetrahedra_element_metrics.py
  // check barycentric coordinates of the point relative to this
  // triangle. If we violate barycentric coorindate assumptions
  // (l1234 >= 0, l1+l2+l3+l4 = 1) we are not inside the hull.
  double l1 = (x*y2*z3 - x*y2*z4 - x*y3*z2 + x*y3*z4 + x*y4*z2 - x*y4*z3 - x2*y*z3 + x2*y*z4 + x2*y3*z - x2*y3*z4 - x2*y4*z + x2*y4*z3 + x3*y*z2 - x3*y*z4 - x3*y2*z + x3*y2*z4 + x3*y4*z - x3*y4*z2 - x4*y*z2 + x4*y*z3 + x4*y2*z - x4*y2*z3 - x4*y3*z + x4*y3*z2)/(x1*y2*z3 - x1*y2*z4 - x1*y3*z2 + x1*y3*z4 + x1*y4*z2 - x1*y4*z3 - x2*y1*z3 + x2*y1*z4 + x2*y3*z1 - x2*y3*z4 - x2*y4*z1 + x2*y4*z3 + x3*y1*z2 - x3*y1*z4 - x3*y2*z1 + x3*y2*z4 + x3*y4*z1 - x3*y4*z2 - x4*y1*z2 + x4*y1*z3 + x4*y2*z1 - x4*y2*z3 - x4*y3*z1 + x4*y3*z2);
  double l2 = (-x*y1*z3 + x*y1*z4 + x*y3*z1 - x*y3*z4 - x*y4*z1 + x*y4*z3 + x1*y*z3 - x1*y*z4 - x1*y3*z + x1*y3*z4 + x1*y4*z - x1*y4*z3 - x3*y*z1 + x3*y*z4 + x3*y1*z - x3*y1*z4 - x3*y4*z + x3*y4*z1 + x4*y*z1 - x4*y*z3 - x4*y1*z + x4*y1*z3 + x4*y3*z - x4*y3*z1)/(x1*y2*z3 - x1*y2*z4 - x1*y3*z2 + x1*y3*z4 + x1*y4*z2 - x1*y4*z3 - x2*y1*z3 + x2*y1*z4 + x2*y3*z1 - x2*y3*z4 - x2*y4*z1 + x2*y4*z3 + x3*y1*z2 - x3*y1*z4 - x3*y2*z1 + x3*y2*z4 + x3*y4*z1 - x3*y4*z2 - x4*y1*z2 + x4*y1*z3 + x4*y2*z1 - x4*y2*z3 - x4*y3*z1 + x4*y3*z2);
  double l3 = ((x - x4)*((y1 - y4)*((x1 - x4)*(-z2 + z4) + (x2 - x4)*(z1 - z4)) + (z1 - z4)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4))) + (x1 - x4)*(-y*((x1 - x4)*(-z2 + z4) + (x2 - x4)*(z1 - z4)) + y4*((x1 - x4)*(-z2 + z4) + (x2 - x4)*(z1 - z4)) - z*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4)) + z4*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4))))/(-(x1 - x4)*(z3 - z4)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4)) + (x3 - x4)*(z1 - z4)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4)) + ((x1 - x4)*(y3 - y4) - (x3 - x4)*(y1 - y4))*((x1 - x4)*(z2 - z4) - (x2 - x4)*(z1 - z4)));
  double l4 = 1-l1-l2-l3;
    
  if(l1 < 0) return false;
  if(l2 < 0) return false;
  if(l3 < 0) return false;
  if(l4 < 0) return false;
    
  if(fabs(1-(l1+l2+l3)) < 1e-8) return false;

  // no assumptions violated, inside triangle
  return true;
    
}

Eigen::Vector4d TetP1::interpolateAtPoint(double r, double s, double t) {
  // barycentric coordinates l1,l2,l3,l4 (l1+l2+l3+l4=1) give us a linear weighting between coordinates (r,s,t). By computing the barycentric weights, we can compute the interpolated velocity = [l1,l2,l3].dot([v1,v2,v3]) where vN is the vertex velocity.
  // T*[l1;l2;l3] = xyz - v4, where xyz is the desired point, and v4 is the vertex 4
  // thus l123 = T**(-1)*xyz - T**(-1)*v4
  // For reference tet
  // Tref =
  // Matrix([
  // [ 0,  0,  1],
  // [ 0,  1,  0],
  // [-1, -1, -1]])
  // See tetrahedra_element_metrics.py: `interp_vector`
  Eigen::Vector4d interpolator;
  interpolator <<
    -r - s - t + 1, s, r, t;    
  return interpolator;
}

Eigen::Vector3d TetP1::inverseCoordinateTransform(double x, double y, double z,
                                                  const Ref<const Matrix<double,mNumVtx,mNumDim>> &vtx) {

  auto x1 = vtx(0,0);
  auto x2 = vtx(1,0);
  auto x3 = vtx(2,0);
  auto x4 = vtx(3,0);
            
  auto y1 = vtx(0,1);
  auto y2 = vtx(1,1);
  auto y3 = vtx(2,1);
  auto y4 = vtx(3,1);
            
  auto z1 = vtx(0,2);
  auto z2 = vtx(1,2);
  auto z3 = vtx(2,2);
  auto z4 = vtx(3,2);

  // see tetrahedra_element_metrics.py (and https://en.wikipedia.org/wiki/Barycentric_coordinate_system)
  auto r=((x - x4)*(-(y1 - y4)*((x1 - x4)*(z2 - z4) - (x2 - x4)*(z1 - z4)) + (z1 - z4)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4))) - (x1 - x4)*(y*(-(x1 - x4)*(z2 - z4) + (x2 - x4)*(z1 - z4)) - y4*(-(x1 - x4)*(z2 - z4) + (x2 - x4)*(z1 - z4)) + z*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4)) - z4*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4))))/(-(x1 - x4)*(z3 - z4)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4)) + (x3 - x4)*(z1 - z4)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4)) + ((x1 - x4)*(y3 - y4) - (x3 - x4)*(y1 - y4))*((x1 - x4)*(z2 - z4) - (x2 - x4)*(z1 - z4)));

  auto s=(-x*y1*z3 + x*y1*z4 + x*y3*z1 - x*y3*z4 - x*y4*z1 + x*y4*z3 + x1*y*z3 - x1*y*z4 - x1*y3*z + x1*y3*z4 + x1*y4*z - x1*y4*z3 - x3*y*z1 + x3*y*z4 + x3*y1*z - x3*y1*z4 - x3*y4*z + x3*y4*z1 + x4*y*z1 - x4*y*z3 - x4*y1*z + x4*y1*z3 + x4*y3*z - x4*y3*z1)/(x1*y2*z3 - x1*y2*z4 - x1*y3*z2 + x1*y3*z4 + x1*y4*z2 - x1*y4*z3 - x2*y1*z3 + x2*y1*z4 + x2*y3*z1 - x2*y3*z4 - x2*y4*z1 + x2*y4*z3 + x3*y1*z2 - x3*y1*z4 - x3*y2*z1 + x3*y2*z4 + x3*y4*z1 - x3*y4*z2 - x4*y1*z2 + x4*y1*z3 + x4*y2*z1 - x4*y2*z3 - x4*y3*z1 + x4*y3*z2);

  auto t=(-x*y1*z2 + x*y1*z3 + x*y2*z1 - x*y2*z3 - x*y3*z1 + x*y3*z2 + x1*y*z2 - x1*y*z3 - x1*y2*z + x1*y2*z3 + x1*y3*z - x1*y3*z2 - x2*y*z1 + x2*y*z3 + x2*y1*z - x2*y1*z3 - x2*y3*z + x2*y3*z1 + x3*y*z1 - x3*y*z2 - x3*y1*z + x3*y1*z2 + x3*y2*z - x3*y2*z1)/(x1*y2*z3 - x1*y2*z4 - x1*y3*z2 + x1*y3*z4 + x1*y4*z2 - x1*y4*z3 - x2*y1*z3 + x2*y1*z4 + x2*y3*z1 - x2*y3*z4 - x2*y4*z1 + x2*y4*z3 + x3*y1*z2 - x3*y1*z4 - x3*y2*z1 + x3*y2*z4 + x3*y4*z1 - x3*y4*z2 - x4*y1*z2 + x4*y1*z3 + x4*y2*z1 - x4*y2*z3 - x4*y3*z1 + x4*y3*z2);    

  Eigen::Vector3d solution {r, s, t};
  return solution;
}

std::tuple<Eigen::Matrix3d,PetscReal> TetP1::inverseJacobian(const Ref<const Matrix<double,mNumVtx,mNumDim>>
                                                             &vtx) {

  auto v1x = vtx(0,0);
  auto v2x = vtx(1,0);
  auto v3x = vtx(2,0);
  auto v4x = vtx(3,0);
             
  auto v1y = vtx(0,1);
  auto v2y = vtx(1,1);
  auto v3y = vtx(2,1);
  auto v4y = vtx(3,1);
             
  auto v1z = vtx(0,2);
  auto v2z = vtx(1,2);
  auto v3z = vtx(2,2);
  auto v4z = vtx(3,2);
  
  // follows correct vertex ordering (see header for diagram)
  auto dxdr = (v3x-v1x);
  auto dxds = (v2x-v1x);
  auto dxdt = (v4x-v1x);
  auto dydr = (v3y-v1y);
  auto dyds = (v2y-v1y);
  auto dydt = (v4y-v1y);
  auto dzdr = (v3z-v1z);
  auto dzds = (v2z-v1z);
  auto dzdt = (v4z-v1z);
  
  Matrix3d J;
  J <<
    dxdr, dydr, dzdr,
    dxds, dyds, dzds,
    dxdt, dydt, dzdt;

  PetscReal detJ = J.determinant();

  // need to test if this analytical inverse is faster than simply J.inverse();
  
  auto  rx = (1 / detJ) * (dyds*dzdt - dydt*dzds);
  auto  sx = (1 / detJ) * (-dydr*dzdt + dydt*dzdr);
  auto  tx = (1 / detJ) * (dydr*dzds - dyds*dzdr);
  
  auto  ry = (1 / detJ) * (-dxds*dzdt + dxdt*dzds);
  auto  sy = (1 / detJ) * ( dxdr*dzdt - dxdt*dzdr);
  auto  ty = (1 / detJ) * (-dxdr*dzds + dxds*dzdr);
  
  auto  rz = (1 / detJ) * (dxds*dydt - dxdt*dyds);
  auto  sz = (1 / detJ) * (-dxdr*dydt + dxdt*dydr);
  auto  tz = (1 / detJ) * (dxdr*dyds - dxds*dydr);

  Matrix3d inverseJacobian;
  inverseJacobian <<
    rx, sx, tx,
    ry, sy, ty,
    rz, sz, tz;

  return std::make_tuple(inverseJacobian, detJ);
  
}

// global x-z points on all nodes
std::tuple<Eigen::VectorXd,
           Eigen::VectorXd,
           Eigen::VectorXd> TetP1::buildNodalPoints(const Ref<const VectorXd>& intCrdR,
                                                    const Ref<const VectorXd>& intCrdS,
                                                    const Ref<const VectorXd>& intCrdT,
                                                    const Ref<const Matrix<double,mNumVtx,mNumDim>>& vtx) {
		
  std::vector<PetscReal> ni(vtx.size());
	
  Eigen::VectorXd nodalPoints_x(intCrdR.size());
  Eigen::VectorXd nodalPoints_y(intCrdR.size());
  Eigen::VectorXd nodalPoints_z(intCrdR.size());

  auto x1 = vtx(0,0);
  auto x2 = vtx(1,0);
  auto x3 = vtx(2,0);
  auto x4 = vtx(3,0);
            
  auto y1 = vtx(0,1);
  auto y2 = vtx(1,1);
  auto y3 = vtx(2,1);
  auto y4 = vtx(3,1);
            
  auto z1 = vtx(0,2);
  auto z2 = vtx(1,2);
  auto z3 = vtx(2,2);
  auto z4 = vtx(3,2);

  // printf("v0(%f,%f,%f),v1(%f,%f,%f),v2(%f,%f,%f),v3(%f,%f,%f)\n",
  // x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4);
  
  for(auto n = 0; n < intCrdR.size(); n++) {
    auto r = intCrdR[n];
    auto s = intCrdS[n];
    auto t = intCrdT[n];
    // map from reference tet point (r,s,t) to this tet (x,y,z)
    auto xn=r*x3 + s*x2 + t*x4 - x1*(r + s + t - 1);
    auto yn=r*y3 + s*y2 + t*y4 - y1*(r + s + t - 1);
    auto zn=r*z3 + s*z2 + t*z4 - z1*(r + s + t - 1);
    
    nodalPoints_x(n) = xn;
    nodalPoints_y(n) = yn;
    nodalPoints_z(n) = zn;
  }    

  return std::make_tuple(nodalPoints_x,nodalPoints_y,nodalPoints_z);
}
