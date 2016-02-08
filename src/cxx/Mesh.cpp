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

void Mesh::registerFieldVector(const std::string &name, const bool &check_out, const bool &check_in) {

    Vec field_vector_local;
    Vec field_vector_global;
    DMCreateLocalVector(mDistributedMesh, &field_vector_local);
    DMCreateGlobalVector(mDistributedMesh, &field_vector_global);

    double zero = 0.0;
    VecSet(field_vector_local, zero);
    PetscObjectSetName((PetscObject) field_vector_global, name.c_str());

    vec_struct registrar;
    registrar.name = name;
    registrar.check_in = check_in;
    registrar.check_out = check_out;
    registrar.field_locals = field_vector_local;
    registrar.field_globals = field_vector_global;
    mFields[name] = registrar;
//
//    mFieldVectorNames.push_back(name);
//    mFieldVectorCheckin.push_back(check_in);
//    mFieldVectorCheckout.push_back(check_out);
//    mFieldVectorLocals.push_back(field_vector_local);
//    mFieldVectorGlobals.push_back(field_vector_global);

}

void Mesh::checkOutFields() {

    for (auto i = 0; i < mFieldVectorNames.size(); i++) {
        if (mFieldVectorCheckout[i]) {
            DMGlobalToLocalBegin(mDistributedMesh, mFieldVectorGlobals[i],
                                 INSERT_VALUES, mFieldVectorLocals[i]);
        }
    }

    for (auto i = 0; i < mFieldVectorNames.size(); i++) {
        if (mFieldVectorCheckout[i]) {
            DMGlobalToLocalEnd(mDistributedMesh, mFieldVectorGlobals[i],
                               INSERT_VALUES, mFieldVectorLocals[i]);
            mCheckedOutFields.push_back(mFieldVectorNames[i]);
        }
    }

}

void Mesh::getFieldOnElement(Eigen::VectorXd &field, Eigen::VectorXi &closure, const std::string &name,
                             const int &element_number) {

    int itr = 0;
    PetscScalar *val = NULL;
    for (auto i = 0; i < mFieldVectorNames.size(); i++) {
        if (mFieldVectorNames[i] == name) {
            DMPlexVecGetClosure(mDistributedMesh, mMeshSection, mFieldVectorLocals[i], element_number, NULL, &val);
            for (auto j = 0; j < closure.size(); j++) { field(closure(j)) = val[j]; }
            DMPlexVecRestoreClosure(mDistributedMesh, mMeshSection, mFieldVectorLocals[i], element_number, NULL, &val);
            break;
        }

    }

}

void Mesh::setFieldOnElement(const Eigen::VectorXd &field, Eigen::VectorXi &closure, const std::string &name,
                             const int &element_number) {

    int itr = 0;
    Eigen::VectorXd val(closure.size());
    for (auto i = 0; i < mFieldVectorNames.size(); i++) {
        if (mFieldVectorNames[i] == name) {
            for (auto j = 0; j < closure.size(); j++) { val(j) = field(closure(j)); }
            DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFieldVectorLocals[i], element_number, val.data(),
                                ADD_VALUES);
            break;
        }
    }

}

void Mesh::checkInFieldsBegin() {

    for (auto i = 0; i < mFieldVectorNames.size(); i++) {
        if (mFieldVectorCheckin[i]) {
            DMLocalToGlobalBegin(mDistributedMesh, mFieldVectorLocals[i], ADD_VALUES, mFieldVectorGlobals[i]);
        }
    }

}


void Mesh::checkInFieldsEnd() {

    for (auto i = 0; i < mFieldVectorNames.size(); i++) {
        if (mFieldVectorCheckin[i]) {
            DMLocalToGlobalEnd(mDistributedMesh, mFieldVectorLocals[i], ADD_VALUES, mFieldVectorGlobals[i]);
        }
    }


}

void ScalarNewmark::advanceField() {

    double dt = 1.0;
    double pre_factor_acceleration = 1.0;
    double pre_factor_displacement = 1.0;

//    VecPointwiseMult(mMassMatrix, mFields["acceleration"].field_globals, mFields["force"].field_globals);
//
//    VecAXPBYPCZ(mFields["velocity"].field_globals, pre_factor_acceleration, pre_factor_acceleration, 1.0,
//                mFields["acceleration"].field_globals, mFields["acceleration_"].field_globals);
//
    VecAXPBYPCZ(mFields["displacement"].field_globals, dt, pre_factor_displacement, 1.0,
                mFields["velocity"].field_globals, mFields["acceleration"].field_globals);
//
//    VecCopy(mFields["acceleration"].field_globals, mFields["acceleration_"].field_globals);

}
