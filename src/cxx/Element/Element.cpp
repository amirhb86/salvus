//
// Created by Michael Afanasiev on 2016-01-29.
//

#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <petscdm.h>
#include <petscdmplex.h>
#include "../Utilities.h"
#include "../Options.h"

#include "Element.h"
#include "HyperCube/Square/Acoustic.h"


Element *Element::factory(Options options) {

    std::string element_shape(options.ElementShape());
    std::string physics_system(options.PhysicsSystem());
    PetscInt polynomial_order = options.PolynomialOrder();
    try {
        if (element_shape == "square") {
            if (physics_system == "acoustic") {
                if (polynomial_order == 4) {
                    return new Acoustic(options);
                }
            }
        }
    throw std::runtime_error("Runtime Error: Element type not supported");
    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI::COMM_WORLD.Abort(-1);
        return nullptr;
    }

}



void Element::registerMesh(DM &distributed_mesh, PetscSection &section) {

    mDistributedMesh = distributed_mesh;
    mMeshSection = section;

}

void Element::printInfoToScreen() const {

    std::cout << "ELEMENT INFORMATION"
            << "\nDimension: " << mNumberDimensions
            << "\nDegrees of freedom per vertex: " << mNumberDofVertex
            << "\nDegrees of freedom per edge: " << mNumberDofEdge
            << "\nDegrees of freedom per face: " << mNumberDofFace
            << "\nDegrees of freedom in volume: " << mNumberDofVolume << std::endl;

//    std::cout << "VERTICES" << std::endl;
//    for (auto &vertex: mVertexCoordinates) { std::cout << vertex << std::endl; }
//
//    std::cout << "INTEGRATION POINTS" << std::endl;
//    for (auto &point: mIntegrationPoints) { std::cout << point << std::endl; }

}

Eigen::VectorXd Element::GllPointsForOrder(const int order) {

    Eigen::VectorXd gll_points;
    if (order == 4) {
        gll_points.resize(5);
        gll_points << -1.0, -0.6546536707, 0.0, 0.6546536707, 1.0;
    }

    return gll_points;

}

Eigen::VectorXd Element::GllIntegrationWeightForOrder(const int order) {

    Eigen::VectorXd integration_weights;
    if (order == 4) {
        integration_weights.resize(5);
        integration_weights << 0.1, 0.5444444444, 0.7111111111, 0.5444444444, 0.1;
    }

    return integration_weights;

}

Eigen::VectorXi Element::ClosureMapping(const int order, const int dimension) {

    Eigen::VectorXi closure_mapping;
    if (dimension == 2) {
        if (order == 4) {
            closure_mapping.resize(25);
            closure_mapping << 6, 13, 22, 3, 15, 7, 16, 23, 2, 20, 8, 17,
                    19, 1, 24, 11, 18, 14, 5, 4, 12, 21, 9, 10, 0;
        }
    }

    return closure_mapping;
}
