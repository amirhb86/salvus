//
// Created by Michael Afanasiev on 2016-01-29.
//

#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <petscdm.h>
#include <petscdmplex.h>
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

void Element::__gatherDistributedFieldsToPartition(Vec &local_field_vec, Vec &global_field_vec) {

    DMGlobalToLocalBegin(mDistributedMesh, global_field_vec, INSERT_VALUES, local_field_vec);
    DMGlobalToLocalEnd(mDistributedMesh, global_field_vec, INSERT_VALUES, local_field_vec);

}

void Element::__gatherPartitionFieldsToElement(Vec &local_field_vec,
                                               Eigen::VectorXd &element_field_vec,
                                               std::vector<int> &closure_mapping) {

    int itr = 0;
    PetscScalar  *val = NULL;
    DMPlexVecGetClosure(mDistributedMesh, mMeshSection, local_field_vec, mLocalElementNumber, NULL, &val);
    for (auto &i: closure_mapping) { element_field_vec(i) = val[itr]; itr++; }
    DMPlexVecRestoreClosure(mDistributedMesh, mMeshSection, local_field_vec, mLocalElementNumber, NULL, &val);

}

void Element::__scatterElementFieldsToPartition(Vec &local_field_vec, const Eigen::VectorXd &element_field_vec) {
    DMPlexVecSetClosure(mDistributedMesh, mMeshSection, local_field_vec, mLocalElementNumber,
                        element_field_vec.data(), ADD_VALUES);
}

void Element::__scatterPartitionFieldsToDistributedBegin(Vec &local_field_vec, Vec &global_field_vec) {
    DMLocalToGlobalBegin(mDistributedMesh, local_field_vec, ADD_VALUES, global_field_vec);
}

void Element::__scatterPartitionFieldsToDistributedEnd(Vec &local_field_vec, Vec &global_field_vec) {
    DMLocalToGlobalEnd(mDistributedMesh, local_field_vec, ADD_VALUES, global_field_vec);
}

