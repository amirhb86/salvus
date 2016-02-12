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
    VecSet(field_vector_global, zero);
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
    double vecmax;
    VecMax(mFields[field_num].field_locals, NULL, &vecmax);

    for (auto j = 0; j < closure.size(); j++) { val(j) = field(closure(j)); }
    DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[field_num].field_locals,
                        element_number, val.data(), ADD_VALUES);

}

void Mesh::checkInFieldsBegin() {

    for (auto reg: mFieldVectorCheckin) {
        double vecmax;
        VecMax(mFields[reg].field_locals, NULL, &vecmax);
        DMLocalToGlobalBegin(mDistributedMesh, mFields[reg].field_locals, ADD_VALUES, mFields[reg].field_globals);
    }

}


void Mesh::checkInFieldsEnd() {

    for (auto reg: mFieldVectorCheckin) {
        DMLocalToGlobalEnd(mDistributedMesh, mFields[reg].field_locals, ADD_VALUES, mFields[reg].field_globals);
    }

}

void Mesh::checkInMassMatrix() {

    int index = -1;
    DMLocalToGlobalBegin(mDistributedMesh, mFields[index].field_locals, ADD_VALUES, mFields[index].field_globals);
    DMLocalToGlobalEnd(mDistributedMesh, mFields[index].field_locals, ADD_VALUES, mFields[index].field_globals);

}

void ScalarNewmark::advanceField() {

    int acceleration_ = (int) AcousticFields::acceleration_;
    int acceleration = (int) AcousticFields::acceleration;
    int displacement = (int) AcousticFields::displacement;
    int velocity = (int) AcousticFields::velocity;

    double dt = 1e-3;
    double pre_factor_acceleration = (1.0/2.0) * dt;
    double pre_factor_displacement = (1.0/2.0) * (dt * dt);

    VecAXPBYPCZ(mFields[velocity].field_globals, pre_factor_acceleration, pre_factor_acceleration, 1.0,
                mFields[acceleration].field_globals, mFields[acceleration_].field_globals);

    VecAXPBYPCZ(mFields[displacement].field_globals, dt, pre_factor_displacement, 1.0,
                mFields[velocity].field_globals, mFields[acceleration].field_globals);

    VecCopy(mFields[acceleration].field_globals, mFields[acceleration_].field_globals);

}

void ScalarNewmark::applyInverseMassMatrix() {

    int mass_matrix_inverse = (int) AcousticFields::mass_matrix_inverse;
    int acceleration = (int) AcousticFields::acceleration;
    int mass_matrix = (int) AcousticFields::mass_matrix;
    int force = (int) AcousticFields::force;

    if (mFields.find(mass_matrix_inverse) == mFields.end()){
        registerFieldVectors(mass_matrix_inverse, false, false, "mass_matrix_inverse");
        VecCopy(mFields[mass_matrix].field_globals, mFields[mass_matrix_inverse].field_globals);
        VecReciprocal(mFields[mass_matrix_inverse].field_globals);
    }

    VecPointwiseMult(mFields[acceleration].field_globals, mFields[mass_matrix_inverse].field_globals,
                     mFields[force].field_globals);

}

void Mesh::zeroFields() {

    double zero = 0.0;
    for (auto reg: mFieldVectorCheckin) {
        VecSet(mFields[reg].field_locals, zero);
        VecSet(mFields[reg].field_globals, zero);
    }

}

void Mesh::setUpMovie() {

    mTime = 0;
    mViewer = nullptr;
    std::string file_name = "/Users/michaelafanasiev/Desktop/movie.h5";
    PetscViewerHDF5Open(PETSC_COMM_WORLD, file_name.c_str(), FILE_MODE_WRITE, &mViewer);
    PetscViewerHDF5PushGroup(mViewer, "/");
    DMView(mDistributedMesh, mViewer);

}

void Mesh::saveFrame() {

    DMSetOutputSequenceNumber(mDistributedMesh, mTime, mTime);
    VecView(mFields[(int) AcousticFields::displacement].field_globals, mViewer);
//    VecView(mFields[-1].field_globals, mViewer);
    mTime += 1;

}

void Mesh::finalizeMovie() {

    PetscViewerHDF5PopGroup(mViewer);
    PetscViewerDestroy(&mViewer);

}
