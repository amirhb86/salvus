//
// Created by Michael Afanasiev on 2016-01-27.
//

#include <mpi.h>
#include <assert.h>
#include "Mesh.h"
#include "Utilities.h"
#include "petscdm.h"
#include "petscdmplex.h"

Mesh *Mesh::factory(Options options) {

    std::string mesh_type(options.MeshType());
    try {
        if (mesh_type == "exodus") {
            return new Exodus;
        } else {
            throw std::runtime_error("Runtime Error: Mesh type + " + mesh_type + " not supported");
        }
    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI::COMM_WORLD.Abort(-1);
        return nullptr;
    };

}

void Mesh::read(Options options) {

    // Class variables.
    mDistributedMesh = NULL;
    mExodusFileName = options.ExodusMeshFile();

    // Function variables.
    DM dm = NULL;
    PetscInt partition_overlap = 0;
    PetscBool interpolate_edges = PETSC_TRUE;

    // Read exodus file.
    DMPlexCreateExodusFromFile(PETSC_COMM_WORLD, mExodusFileName.c_str(), interpolate_edges, &dm);

    // May be a race condition on distribute.
    MPI::COMM_WORLD.Barrier();
    DMPlexDistribute(dm, partition_overlap, NULL, &mDistributedMesh);

    // We don't need the serial mesh anymore.
    if (mDistributedMesh) { DMDestroy(&dm); }

    // Set some important information.
    DMGetDimension(mDistributedMesh, &mNumberDimensions);
    DMPlexGetDepthStratum(mDistributedMesh, mNumberDimensions, NULL, &mNumberElementsLocal);

}

void Mesh::setupGlobalDof(PetscInt number_dof_vertex, PetscInt number_dof_edge, PetscInt number_dof_face,
                          PetscInt number_dof_volume, PetscInt number_dimensions) {

    // Ensure that the mesh and the elements are at least the same dimension.
    assert(number_dimensions == mNumberDimensions);

    // Only define 1 field here because we're taking care of multiple fields manually.
    PetscInt number_fields = 1;
    PetscInt number_components = 1;
    PetscInt number_dof_per_element[mNumberDimensions + 1];

    // Number of dof on vertex, edge, face, volume.
    number_dof_per_element[0] = number_dof_vertex;
    number_dof_per_element[1] = number_dof_edge;
    number_dof_per_element[2] = number_dof_face;
    if (mNumberDimensions == 3) { number_dof_per_element[3] = number_dof_volume; }

    // Setup the global and local (distributed) degrees of freedom.
    DMPlexCreateSection(mDistributedMesh, mNumberDimensions, number_fields, &number_components,
                        number_dof_per_element, 0, NULL, NULL, NULL, NULL, &mMeshSection);
    DMSetDefaultSection(mDistributedMesh, mMeshSection);

}

std::vector<PetscReal> Mesh::ElementVertices(const PetscInt element_number) {


}
