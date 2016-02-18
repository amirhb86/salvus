//
// Created by Michael Afanasiev on 2016-01-27.
//

#include <mpi.h>
#include <petscviewerhdf5.h>
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
            return new ScalarNewmark2D;
        } else {
            throw std::runtime_error("Runtime Error: Mesh type " + mesh_type + " not supported");
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

void Mesh::registerFieldVectors(const std::string &name) {

    Vec field_vector_local;
    Vec field_vector_global;
    DMCreateLocalVector(mDistributedMesh, &field_vector_local);
    DMCreateGlobalVector(mDistributedMesh, &field_vector_global);

    double zero = 0.0;
    VecSet(field_vector_local, zero);
    VecSet(field_vector_global, zero);
    PetscObjectSetName((PetscObject) field_vector_global, name.c_str());

    vec_struct registrar;
    registrar.name = name;
    registrar.loc = field_vector_local;
    registrar.glb = field_vector_global;
    mFields[name] = registrar;

}

void Mesh::checkOutField(const std::string &name) {

    assert(mFields.find(name) != mFields.end());

    // Begin the MPI broadcast global -> local.
    DMGlobalToLocalBegin(mDistributedMesh, mFields[name].glb, INSERT_VALUES,
                         mFields[name].loc);

    // End the MPI broadcast global -> local.
    DMGlobalToLocalEnd(mDistributedMesh, mFields[name].glb, INSERT_VALUES,
                       mFields[name].loc);

}

Eigen::VectorXd Mesh::getFieldOnElement(const std::string &name, const int &element_number,
                                        const Eigen::VectorXi &closure) {

    PetscScalar *val = NULL;
    Eigen::VectorXd field(closure.size());
    DMPlexVecGetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                        element_number, NULL, &val);
    for (auto j = 0; j < closure.size(); j++) { field(closure(j)) = val[j]; }
    DMPlexVecRestoreClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                            element_number, NULL, &val);
    return field;

}

void Mesh::setFieldOnElement(const std::string &name, const int &element_number,
                             const Eigen::VectorXi &closure, const Eigen::VectorXd &field) {

    Eigen::VectorXd val(closure.size());
    // map "our" nodal ordering back onto PETSC ordering
    for (auto j = 0; j < closure.size(); j++) { val(j) = field(closure(j)); }
    DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                        element_number, val.data(), ADD_VALUES);

}

void Mesh::setFieldFromElement(const std::string &name, const int element_number,
                               const Eigen::VectorXi &closure, const Eigen::VectorXd &field) {

    Eigen::VectorXd val(closure.size());
    // map "our" nodal ordering back onto PETSC ordering
    for (auto j = 0; j < closure.size(); j++) { val(j) = field(closure(j)); }
    DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                        element_number, val.data(), ADD_VALUES);
}

void Mesh::addFieldFromElement(const std::string &name, const int element_number,
                               const Eigen::VectorXi &closure, const Eigen::VectorXd &field) {

    Eigen::VectorXd val(closure.size());
    // map "our" nodal ordering back onto PETSC ordering
    for (auto j = 0; j < closure.size(); j++) { val(j) = field(closure(j)); }
    DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                        element_number, val.data(), ADD_VALUES);
}

void Mesh::assembleLocalFieldToGlobal(const std::string &name) {

    assembleLocalFieldToGlobalBegin(name);
    assembleLocalFieldToGlobalEnd(name);   
}

void Mesh::assembleLocalFieldToGlobalBegin(const std::string &name) {

    // Make sure the field exists in our dictionary.
    assert(mFields.find(name) != mFields.end());

    // Begin MPI broadcast local -> global.
    DMLocalToGlobalBegin(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);
}

void Mesh::assembleLocalFieldToGlobalEnd(const std::string &name) {

    // Make sure the field exists in our dictionary.
    assert(mFields.find(name) != mFields.end());

    // Begin MPI broadcast local -> global.
    DMLocalToGlobalEnd(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);
}

void Mesh::setLocalFieldToGlobal(const std::string &name) {
    
    // Make sure the field exists in our dictionary.
    assert(mFields.find(name) != mFields.end());

    // Do "communication". `INSERT_VALUE` will result in no communication
//    DMLocalToGlobalBegin(mDistributedMesh, mFields[name].loc, INSERT_VALUES, mFields[name].glb);
//    DMLocalToGlobalEnd(mDistributedMesh, mFields[name].loc, INSERT_VALUES, mFields[name].glb);
    DMLocalToGlobalBegin(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);
    DMLocalToGlobalEnd(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);

}

// Depricated
void Mesh::checkInFieldBegin(const std::string &name) {

    std::cout << "WARNING: `checkInFieldBegin` To be depricated\n";
    // Make sure the field exists in our dictionary.
    assert(mFields.find(name) != mFields.end());

    // Begin MPI broadcast local -> global.
    DMLocalToGlobalBegin(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);
}


void Mesh::checkInFieldEnd(const std::string &name) {

    std::cout << "WARNING: `checkInFieldEnd` To be depricated\n";
    // Make sure the field exists in our dictionary.
    assert(mFields.find(name) != mFields.end());

    // Begin MPI broadcast local -> global.
    DMLocalToGlobalEnd(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);
}

void Mesh::zeroFields(const std::string &name) {
    double zero = 0.0;
    VecSet(mFields[name].loc, zero);
    VecSet(mFields[name].glb, zero);
}

void Mesh::setUpMovie(const std::string &movie_filename) {
    mTime = 0;
    mViewer = nullptr;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, movie_filename.c_str(), FILE_MODE_WRITE, &mViewer);
    PetscViewerHDF5PushGroup(mViewer, "/");
    DMView(mDistributedMesh, mViewer);
}

void Mesh::saveFrame(std::string name) {

    DMSetOutputSequenceNumber(mDistributedMesh, mTime, mTime);
    VecView(mFields[name].glb, mViewer);
    mTime += 1;

}

void Mesh::finalizeMovie() {

    PetscViewerHDF5PopGroup(mViewer);
    PetscViewerDestroy(&mViewer);

}

void ScalarNewmark2D::advanceField() {

    double dt = 1e-3;
    double pre_factor_acceleration = (1.0/2.0) * dt;
    double pre_factor_displacement = (1.0/2.0) * (dt * dt);

    VecAXPBYPCZ(mFields["velocity"].glb, pre_factor_acceleration, pre_factor_acceleration, 1.0,
                mFields["acceleration"].glb, mFields["acceleration_"].glb);

    VecAXPBYPCZ(mFields["displacement"].glb, dt, pre_factor_displacement, 1.0,
                mFields["velocity"].glb, mFields["acceleration"].glb);

    VecCopy(mFields["acceleration"].glb, mFields["acceleration_"].glb);

}

void ScalarNewmark2D::applyInverseMassMatrix() {

    if (mFields.find("mass_matrix_inverse") == mFields.end()){
        registerFieldVectors("mass_matrix_inverse");
        VecCopy(mFields["mass_matrix"].glb, mFields["mass_matrix_inverse"].glb);
        VecReciprocal(mFields["mass_matrix_inverse"].glb);
    }

    double maxval;
    VecMax(mFields["force"].glb, NULL, &maxval);
    std::cout << "MAX FORCE " << maxval << std::endl;
    VecPointwiseMult(mFields["acceleration"].glb, mFields["mass_matrix_inverse"].glb,
                     mFields["force"].glb);

}

void ElasticNewmark2D::advanceField() {

    double dt = 1e-3;
    double pre_factor_acceleration = (1.0/2.0) * dt;
    double pre_factor_displacement = (1.0/2.0) * (dt * dt);

    VecAXPBYPCZ(mFields["velocity_x"].glb, pre_factor_acceleration, pre_factor_acceleration, 1.0,
                mFields["acceleration_x"].glb, mFields["acceleration_x_"].glb);
    VecAXPBYPCZ(mFields["velocity_z"].glb, pre_factor_acceleration, pre_factor_acceleration, 1.0,
                mFields["acceleration_z"].glb, mFields["acceleration_z_"].glb);

    VecAXPBYPCZ(mFields["displacement_x"].glb, dt, pre_factor_displacement, 1.0,
                mFields["velocity_x"].glb, mFields["acceleration_x"].glb);
    VecAXPBYPCZ(mFields["displacement_z"].glb, dt, pre_factor_displacement, 1.0,
                mFields["velocity_z"].glb, mFields["acceleration_z"].glb);

    VecCopy(mFields["acceleration_x"].glb, mFields["acceleration_x_"].glb);
    VecCopy(mFields["acceleration_z"].glb, mFields["acceleration_z_"].glb);

}

void ElasticNewmark2D::applyInverseMassMatrix() {

    if (mFields.find("mass_matrix_inverse") == mFields.end()) {
        registerFieldVectors("mass_matrix_inverse");
        VecCopy(mFields["mass_matrix"].glb, mFields["mass_matrix_inverse"].glb);
        VecReciprocal(mFields["mass_matrix_inverse"].glb);
    }

    VecPointwiseMult(mFields["acceleration_x"].glb, mFields["mass_matrix_inverse"].glb,
                     mFields["force_x"].glb);
    VecPointwiseMult(mFields["acceleration_z"].glb, mFields["mass_matrix_inverse"].glb,
                     mFields["force_z"].glb);

}
