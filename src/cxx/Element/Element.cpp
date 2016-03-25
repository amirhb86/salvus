#include "Element.h"

#include <Eigen/Dense>

#include "Element/HyperCube/Quad/Acoustic.h"
#include "Element/HyperCube/Quad/Elastic.h"

#include "Element/Simplex/Triangle/AcousticTri.h"

/*
 * STATIC variables WHICH ARE ONLY ON THE REFERENCE ELEMENT.
 */
int Element2D::mNumberDofVertex;
int Element2D::mNumberDofEdge;
int Element2D::mNumberDofFace;
// -1 as a canary --- the element constants haven't been assigned yet
int Element2D::mNumberIntegrationPoints = -1; 
int Element2D::mPolynomialOrder;
Eigen::VectorXi Element2D::mClosureMapping;


Element2D *Element2D::factory(Options options) {

    std::string physics(options.PhysicsSystem());
    
    try {
        if(options.ElementShape() == "quad") {
            std::cout << "Building quad!\n";
            if (physics == "acoustic") {
                return new AcousticQuad(options);
            } else if (physics == "elastic") {
                return new Elastic(options);
            } else {
                throw std::runtime_error("Runtime Error: Element physics " + physics + " not supported");
            }
        }
        else if(options.ElementShape() == "triangle") {
            if (physics == "acoustic") {
                std::cout << "New AcousticTri\n";
                return new AcousticTri(options);
            }
            else {
                PRINT_ROOT() << "ERROR: Physics for triangles not implemented\n";
                MPI::COMM_WORLD.Abort(-1);
                return nullptr;
            }
        }
    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI::COMM_WORLD.Abort(-1);
        return nullptr;
    }
}

Eigen::VectorXd Element2D::checkOutFieldElement(Mesh *mesh, const std::string name) {

    return mesh->getFieldOnElement(name, mElementNumber, mClosureMapping);

}

void Element2D::checkInFieldElement(Mesh *mesh, const Eigen::VectorXd &field, const std::string name) {

    mesh->addFieldFromElement(name, mElementNumber, mClosureMapping, field);

}

void Element2D::setFieldElement(Mesh *mesh, const Eigen::VectorXd &field, const std::string name) {

    mesh->setFieldFromElement(name, mElementNumber, mClosureMapping, field);

}

void Element2D::setBoundaryConditions(Mesh *mesh) {

    mOnBoundary = false;
    for(auto& keys : mesh->BoundaryElementFaces()) {
        auto boundary_name = keys.first;
        auto elements_in_boundary = keys.second;
        if(elements_in_boundary.find(mElementNumber) != elements_in_boundary.end()) {
            // found this element on some boundary            
            mOnBoundary = true;
            // assign boundary_name -> {face_ids}
            mBoundaries[boundary_name] = elements_in_boundary[mElementNumber];
        }
    }
}

void Element2D::applyBoundaryConditions(Mesh *mesh,
                                   Options &options,
                                   std::string fieldname) {

    if(!mOnBoundary) return;
    
    // dirichlet boundaries
    double value = 0; // value to set
    auto dirichlet_boundary_names = options.DirichletBoundaries();
    for(auto& bndry : dirichlet_boundary_names) {
        auto faceids = mBoundaries[bndry];
        for(auto& faceid : faceids) {
            auto field = mesh->getFieldOnFace(fieldname,faceid);
            // apply dirichlet condition
            printf("Setting 0\n");
            field = 0*field.array() + value;
            mesh->setFieldFromFace(fieldname,faceid,field);
        }
    }
}
