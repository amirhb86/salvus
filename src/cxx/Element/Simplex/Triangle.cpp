
#include <mpi.h>
#include <iostream>
#include <petscdm.h>
#include <petscdmplex.h>
#include <tuple>
#include "Triangle.h"

/*
 * STATIC variables WHICH ARE ONLY ON THE REFERENCE ELEMENT.
 */
Eigen::VectorXi Triangle::mClosureMapping;
Eigen::MatrixXd Triangle::mGradientOperator;
Eigen::VectorXd Triangle::mIntegrationWeights;
Eigen::VectorXd Triangle::mIntegrationCoordinates;

Triangle::Triangle(Options options) {

    // mVertexCoordinates has 3 vertices
    mVertexCoordinates.resize(2,mNumberVertex);
    
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






