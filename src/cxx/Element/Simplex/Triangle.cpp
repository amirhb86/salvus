#include <iostream>
#include <petscdm.h>
#include <petscdmplex.h>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <tuple>
#include "Triangle.h"


/*
 * STATIC variables WHICH ARE ONLY ON THE REFERENCE ELEMENT.
 */

Eigen::MatrixXd Triangle::mGradientOperator;
Eigen::VectorXd Triangle::mIntegrationWeights;
Eigen::VectorXd Triangle::mIntegrationCoordinates_r;
Eigen::VectorXd Triangle::mIntegrationCoordinates_s;
Eigen::MatrixXd Triangle::mGradientPhi_dr;
Eigen::MatrixXd Triangle::mGradientPhi_ds;

    
std::tuple<Eigen::VectorXd,Eigen::VectorXd> Triangle::QuadraturePointsForOrder(const int order) {

    if(order == 3) {
        int num_pts = mNumberIntegrationPoints;
        Eigen::VectorXd rn(num_pts);
        Eigen::VectorXd sn(num_pts);
        coordinates_p3_triangle_rn(rn.data());
        coordinates_p3_triangle_sn(sn.data());
        return std::make_tuple(rn,sn);
    }
    else {
        std::cerr << "ERROR: Order NOT implemented!...\n";
        MPI::COMM_WORLD.Abort(-1);
    }
    
}

Eigen::VectorXd Triangle::QuadratureIntegrationWeightForOrder(const int order) {
    
    if(order == 3) {
        int num_pts = mNumberIntegrationPoints;
        Eigen::VectorXd wn(num_pts);
        quadrature_weights_p3_triangle(wn.data());
        return wn;
    } else {
        std::cerr << "ERROR: Order NOT implemented!\n";
        MPI::COMM_WORLD.Abort(-1);
    }
}

Eigen::VectorXi Triangle::ClosureMapping(const int order, const int dimension) {

    if(order == 3) {

        int num_pts = mNumberIntegrationPoints;

        Eigen::VectorXi closure(num_pts);
        closure << 0,1,2,3,4,5,6,7,8,9,10,11;
        return closure;
        
    } else {
        std::cerr << "ERROR: Order NOT implemented!\n";
        MPI::COMM_WORLD.Abort(-1);
    }
    
}



Triangle::Triangle(Options options) {

    std::cout << "Options constructor \n";
    // set to -1 at compile time, and now set only once
    
        
        // Basic properties.
        mPolynomialOrder = options.PolynomialOrder();
        if(mPolynomialOrder == 3) {
            mNumberIntegrationPoints = 12;
            mNumberDofEdge = 2;
            mNumberDofFace = 3;
        }
        else {
            std::cerr << "ERROR: Order NOT implemented!\n";
            MPI::COMM_WORLD.Abort(-1);
        }
        
        // mVertexCoordinates has 3 vertices
        mVertexCoordinates.resize(2,mNumberVertex);
    
        // Nodal collocation points on edge and surface
        mNumberDofVertex = 1;
        
        // Integration points and weights
        std::tie(mIntegrationCoordinates_r,mIntegrationCoordinates_s) =
            Triangle::QuadraturePointsForOrder(options.PolynomialOrder());        
        mIntegrationWeights = Triangle::QuadratureIntegrationWeightForOrder(options.PolynomialOrder());
        
        mClosureMapping = Triangle::ClosureMapping(options.PolynomialOrder(), mNumberDimensions);
        setupGradientOperator();
    
}

void Triangle::attachVertexCoordinates(DM &distributed_mesh) {

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

}


std::tuple<Eigen::Matrix2d,PetscReal> Triangle::inverseJacobianAtPoint(PetscReal eps, PetscReal eta) {

    // Current triangles are affine, so eps and eta are ignored as the
    // tranform from reference element is linear.
    
    double v1x = mVertexCoordinates(0,0);
    double v2x = mVertexCoordinates(0,1);
    double v3x = mVertexCoordinates(0,2);
    
    double v1z = mVertexCoordinates(1,0);
    double v2z = mVertexCoordinates(1,1);
    double v3z = mVertexCoordinates(1,2);

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
    
    Eigen::Matrix<double,2,2> J;
    J << dxdr,dzdr
        ,dxds,dzds;

    auto J_inv = 1/J.determinant();
    auto detJ = 1/J_inv;
    
    auto rx = J_inv*dzds;
    auto rz = -J_inv*dxds;
    auto sx = -J_inv*dzdr;
    auto sz = J_inv*dxdr;
    
    Eigen::Matrix<double,2,2> inverseJacobian;
    inverseJacobian << rx, sx,
                       rz, sz;
    return std::make_tuple(inverseJacobian, detJ);
    
}

Eigen::Vector3d Triangle::__interpolateMaterialProperties(ExodusModel *model, std::string parameter_name) {

    Eigen::Vector3d material_at_vertices(mNumberVertex);

    for (auto i = 0; i < mNumberVertex; i++) {
        material_at_vertices(i) = model->getElementalMaterialParameterAtVertex(
                                                                               mElementCenter, parameter_name, i);
    }
    return material_at_vertices;

}

void Triangle::attachSource(std::vector<Source*> sources) {

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

bool Triangle::mCheckHull(double x, double z) {

    double x1 = mVertexCoordinates.row(0)[0];
    double x2 = mVertexCoordinates.row(0)[1];
    double x3 = mVertexCoordinates.row(0)[2];
    double x4 = mVertexCoordinates.row(0)[3];
    double z1 = mVertexCoordinates.row(1)[0];
    double z2 = mVertexCoordinates.row(1)[1];
    double z3 = mVertexCoordinates.row(1)[2];
    double z4 = mVertexCoordinates.row(1)[3];

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

Eigen::Vector2d Triangle::inverseCoordinateTransform(const double &x, const double &z) {

    double x1 = mVertexCoordinates.row(0)[0];
    double x2 = mVertexCoordinates.row(0)[1];
    double x3 = mVertexCoordinates.row(0)[2];
    double z1 = mVertexCoordinates.row(1)[0];
    double z2 = mVertexCoordinates.row(1)[1];
    double z3 = mVertexCoordinates.row(1)[2];

    // see triangle_element_metrics.py (and https://en.wikipedia.org/wiki/Barycentric_coordinate_system)
    
    double r  = ((-x + x3)*(z2 - z3)*((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)) + (x2 - x3)*(z - z3)*((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)) + (x*(z1 - z3) - x3*(z1 - z3) - z*(x1 - x3) + z3*(x1 - x3))*(x*(z1 - z3) - x3*(z1 - z3) - z*(x1 - x3) + z*(x2 - x3) + z3*(x1 - x3) - z3*(x2 - x3) + (-x + x3)*(z2 - z3) + (x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)))/pow((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3),2.0);
    double s = ((-x + x3)*(z2 - z3)*((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)) + (x2 - x3)*(z - z3)*((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)) + (x*(z1 - z3) - x3*(z1 - z3) - z*(x1 - x3) + z3*(x1 - x3))*(x*(z1 - z3) - x3*(z1 - z3) - z*(x1 - x3) + z*(x2 - x3) + z3*(x1 - x3) - z3*(x2 - x3) + (-x + x3)*(z2 - z3) + (x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3)))/pow((x1 - x3)*(z2 - z3) - (x2 - x3)*(z1 - z3),2.0);

    Eigen::Vector2d solution {r, s};
    return solution;
}

void Triangle::setupGradientOperator() {

    if(mPolynomialOrder == 3) {
        mGradientPhi_dr.resize(mNumberIntegrationPoints,mNumberIntegrationPoints);
        mGradientPhi_ds.resize(mNumberIntegrationPoints,mNumberIntegrationPoints);
        dphi_dr_rsn_p3_triangle(mGradientPhi_dr.data());
        dphi_ds_rsn_p3_triangle(mGradientPhi_ds.data());
    } else {        
        std::cerr << "NOT implemented yet!\n";
        MPI::COMM_WORLD.Abort(-1);
    }
}

// global x-z points on all nodes
std::tuple<Eigen::VectorXd,Eigen::VectorXd> Triangle::buildNodalPoints(Mesh* mesh) {
		
    std::vector<PetscReal> ni(mVertexCoordinates.size());
	
  Eigen::VectorXd nodalPoints_x(mNumberIntegrationPoints);
  Eigen::VectorXd nodalPoints_z(mNumberIntegrationPoints);

  double x1 = mVertexCoordinates.row(0)[0];
  double x2 = mVertexCoordinates.row(0)[1];
  double x3 = mVertexCoordinates.row(0)[2];
  double z1 = mVertexCoordinates.row(1)[0];
  double z2 = mVertexCoordinates.row(1)[1];
  double z3 = mVertexCoordinates.row(1)[2];

  for(auto n = 0; n < mNumberIntegrationPoints; n++) {
      auto r = mIntegrationCoordinates_r[n];
      auto s = mIntegrationCoordinates_s[n];
      // map from reference triangle point (r,s) to this triangle (x,z)
      auto xn=-x1*(r + s)/2 + x2*(r + 1)/2 + x3*(s + 1)/2;
      auto zn=-z1*(r + s)/2 + z2*(r + 1)/2 + z3*(s + 1)/2;
      nodalPoints_x(n) = xn;
      nodalPoints_z(n) = zn;
  }
  
  // push nodal locations to shared dofs
  mesh->setFieldFromElement("x", mElementNumber, mClosureMapping, nodalPoints_x);
  mesh->setFieldFromElement("z", mElementNumber, mClosureMapping, nodalPoints_z);

  return std::make_tuple(nodalPoints_x,nodalPoints_z);
}

double Triangle::integrateField(const Eigen::VectorXd &field) {

    double val = 0;
    Eigen::Matrix<double,2,2> inverse_Jacobian;
    double detJ;
    std::tie(inverse_Jacobian,detJ) = inverseJacobianAtPoint(0,0);
    val = detJ*field.dot(mIntegrationWeights);
    return val;

}
