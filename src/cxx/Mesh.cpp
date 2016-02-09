//
// Created by Michael Afanasiev on 2016-01-27.
//

#include <mpi.h>
#include <assert.h>
#include <petscao.h>
#include "Mesh.h"
#include "Utilities.h"
#include "petscdm.h"
#include "petscdmplex.h"

Mesh *Mesh::factory(Options options) {

    std::string mesh_type(options.MeshType());
    try {
        if (mesh_type == "newmark") {
            return new ScalarNewmark;
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

void Mesh::setupGlobalDof(int number_dof_vertex, int number_dof_edge, int number_dof_face,
                          int number_dof_volume, int number_dimensions) {

    // Ensure that the mesh and the elements are at least the same dimension.
    assert(number_dimensions == mNumberDimensions);

    // Only define 1 field here because we're taking care of multiple fields manually.
    int number_fields = 1;
    int number_components = 1;
    int number_dof_per_element[mNumberDimensions + 1];

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

void Mesh::registerFieldVectors(const int &num, const bool &check_out, const bool &check_in,
                                const std::string &name) {

    Vec field_vector_local;
    Vec field_vector_global;
    DMCreateLocalVector(mDistributedMesh, &field_vector_local);
    DMCreateGlobalVector(mDistributedMesh, &field_vector_global);

    double zero = 0.0;
    VecSet(field_vector_local, zero);
    PetscObjectSetName((PetscObject) field_vector_global, name.c_str());

    if (check_in) { mFieldVectorCheckin.push_back(num); }
    if (check_out) { mFieldVectorCheckout.push_back(num); }

    vec_struct registrar;
    registrar.name = name;
    registrar.check_in = check_in;
    registrar.check_out = check_out;
    registrar.field_locals = field_vector_local;
    registrar.field_globals = field_vector_global;
    mFields[num] = registrar;

}

void Mesh::checkOutFields() {

    for (auto &reg: mFieldVectorCheckout) {
        DMGlobalToLocalBegin(mDistributedMesh, mFields[reg].field_globals, INSERT_VALUES,
                             mFields[reg].field_locals);
    }

    for (auto &reg: mFieldVectorCheckout) {
        DMGlobalToLocalEnd(mDistributedMesh, mFields[reg].field_globals, INSERT_VALUES,
                           mFields[reg].field_locals);
    }

}

Eigen::VectorXd Mesh::getFieldOnElement(const int &field_num, const int &element_number,
                                        const Eigen::VectorXi &closure) {

    PetscScalar *val = NULL;
    Eigen::VectorXd field(closure.size());
    DMPlexVecGetClosure(mDistributedMesh, mMeshSection, mFields[field_num].field_locals,
                        element_number, NULL, &val);
    for (auto j = 0; j < closure.size(); j++) { field(closure(j)) = val[j]; }
    DMPlexVecRestoreClosure(mDistributedMesh, mMeshSection, mFields[field_num].field_locals,
                            element_number, NULL, &val);
    return field;

}

void Mesh::setFieldOnElement(const int &field_num, const int &element_number,
                             const Eigen::VectorXi &closure, const Eigen::VectorXd &field) {

    Eigen::VectorXd val(closure.size());
    for (auto j = 0; j < closure.size(); j++) { val(j) = field(closure(j)); }
    DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[field_num].field_locals,
                        element_number, val.data(), ADD_VALUES);

}

void Mesh::checkInFieldsBegin() {

    for (auto reg: mFieldVectorCheckin) {
        DMLocalToGlobalBegin(mDistributedMesh, mFields[reg].field_locals, ADD_VALUES, mFields[reg].field_globals);
    }

}


void Mesh::checkInFieldsEnd() {

    for (auto reg: mFieldVectorCheckin) {
        DMLocalToGlobalEnd(mDistributedMesh, mFields[reg].field_locals, ADD_VALUES, mFields[reg].field_globals);
    }

}

void ScalarNewmark::advanceField() {

    int acceleration_ = (int) AcousticFields::acceleration_;
    int acceleration = (int) AcousticFields::acceleration;
    int displacement = (int) AcousticFields::displacement;
    int velocity = (int) AcousticFields::velocity;
    int force = (int) AcousticFields::force;

    double dt = 1.0;
    double pre_factor_acceleration = 1.0;
    double pre_factor_displacement = 1.0;

//    VecPointwiseMult(mMassMatrix, mFields[acceleration].field_globals,
//                     mFields[force].field_globals);

    VecAXPBYPCZ(mFields[velocity].field_globals, pre_factor_acceleration, pre_factor_acceleration, 1.0,
                mFields[acceleration].field_globals, mFields[acceleration_].field_globals);

    VecAXPBYPCZ(mFields[displacement].field_globals, dt, pre_factor_displacement, 1.0,
                mFields[velocity].field_globals, mFields[acceleration].field_globals);

    VecCopy(mFields[acceleration].field_globals, mFields[acceleration_].field_globals);

}
