//
// Created by Michael Afanasiev on 2016-01-29.
//

#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <petscdm.h>
#include "../Utilities.h"
#include "../Options.h"

#include "Element.h"
#include "HyperCube/Square/SquareAcousticOrderFour.h"


Element *Element::factory(Options options) {

    std::string element_shape(options.ElementShape());
    std::string physics_system(options.PhysicsSystem());
    PetscInt polynomial_order = options.PolynomialOrder();
    try {
        if (element_shape == "square") {
            if (physics_system == "acoustic") {
                if (polynomial_order == 4) {
                    return new SquareAcousticOrderFour(options);
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
