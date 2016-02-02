//
// Created by Michael Afanasiev on 2016-02-01.
//

#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <assert.h>
#include "exodusII.h"
#include "ExodusModel.h"
#include "../Utilities.h"

ExodusModel::ExodusModel(Options options) {

    mExodusFileName = options.ExodusModelFile();

}

void ExodusModel::initializeParallel() {

    // Do all reading from exodus file on rank 0.
    if (!MPI::COMM_WORLD.Get_rank()) {
        getInitialization();
        readCoordinates();
        readNodalVariables();
    }

    // Broadcast all scalars.
    mNumberDimension = utilities::broadcastInt(mNumberDimension);
    mNumberVertices = utilities::broadcastInt(mNumberVertices);
    mNumberElements = utilities::broadcastInt(mNumberElements);
    mNumberElementBlocks = utilities::broadcastInt(mNumberElementBlocks);
    mNumberNodeSets = utilities::broadcastInt(mNumberNodeSets);
    mNumberSideSets = utilities::broadcastInt(mNumberSideSets);

    // Broadcast all vectors.
    mNodalX = utilities::broadcastStdVecFromRoot(mNodalX);
    mNodalY = utilities::broadcastStdVecFromRoot(mNodalY);
    mNodalVariables = utilities::broadcastStdVecFromRoot(mNodalVariables);
    mNodalVariableNames = utilities::broadcastStringVecFromFroot(mNodalVariableNames);

    // Broadcast dimension specific components.
    if (mNumberDimension > 2) { mNodalZ = utilities::broadcastStdVecFromRoot(mNodalZ); }

    // Create KdTree on each processor.
    createKdTree();

}

void ExodusModel::exodusError(const int retval, std::string func_name) {

    if (retval) {
        throw std::runtime_error("Error in exodus function: " + func_name);
    }

}

void ExodusModel::readCoordinates() {

    try {
        mNodalX.resize(mNumberVertices);
        mNodalY.resize(mNumberVertices);
        if (mNumberDimension > 2) { mNodalZ.resize(mNumberVertices); }
        if (mNumberDimension == 2) {
            exodusError(ex_get_coord(
                                mExodusId, mNodalX.data(), mNodalY.data(), NULL),
                        "ex_get_coord");
        } else {
            exodusError(ex_get_coord(
                                mExodusId, mNodalX.data(), mNodalY.data(), mNodalZ.data()),
                        "ex_get_coord");
        }
    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI_Abort(PETSC_COMM_WORLD, -1);
    }

}

void ExodusModel::getInitialization() {

    int io_ws = 0;
    int comp_ws = 8;
    try {
        mExodusId = ex_open(mExodusFileName.c_str(), EX_READ, &comp_ws, &io_ws, &mExodusVersion);
        if (mExodusId < 0) { throw std::runtime_error("Error opening exodus model file."); }

        exodusError(ex_get_init(
                            mExodusId, mExodusTitle, &mNumberDimension, &mNumberVertices, &mNumberElements,
                            &mNumberElementBlocks, &mNumberNodeSets, &mNumberSideSets),
                    "ex_get_init");

    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI_Abort(PETSC_COMM_WORLD, -1);
    }

}

void ExodusModel::createKdTree() {

    // Need to keep the index arrays allocated.
    mKdTreeData.resize(mNodalX.size());

    // Create.
    if (mNumberDimension == 2) {
        mKdTree = kd_create(2);
        for (auto i = 0; i < mNumberVertices; i++) {
            mKdTreeData[i] = i;
            kd_insert(mKdTree, std::vector<double> {mNodalX[i], mNodalY[i]}.data(), &mKdTreeData[i]);
        }
    } else {
        mKdTree = kd_create(3);
        for (auto i = 0; i < mNumberVertices; i++) {
            mKdTreeData[i] = i;
            kd_insert(mKdTree, std::vector<double> {mNodalX[i], mNodalY[i], mNodalZ[i]}.data(), &mKdTreeData[i]);
        }
    }

}

void ExodusModel::readNodalVariables() {

    // Get variables names.
    exodusError(ex_get_var_param(
            mExodusId, "n", &mNumberNodalVariables),
                "ex_get_var_param");
    char *nm[mNumberNodalVariables];
    for (auto i = 0; i < mNumberNodalVariables; i++) { nm[i] = (char *) calloc((MAX_STR_LENGTH+1), sizeof(char)); }
    exodusError(ex_get_var_names(
            mExodusId, "N", mNumberNodalVariables, nm),
                "ex_get_var_names");
    for (auto i = 0; i < mNumberNodalVariables; i++) { mNodalVariableNames.push_back(std::string(nm[i])); }

    // Get variable values.
    int time_step = 1;
    std::vector<double> buffer(mNumberVertices);
    for (auto i = 0; i < mNumberNodalVariables; i++) {
        exodusError(ex_get_nodal_var(
                            mExodusId, time_step, (i+1), mNumberVertices, buffer.data()),
                    "ex_get_nodal_var");
        mNodalVariables.insert(mNodalVariables.begin(), buffer.begin(), buffer.end());
    }
}

PetscReal ExodusModel::getMaterialParameterAtPoint(const std::vector<double> point,
                                                   const std::string parameter_name) {

    // Ensure dimensions are consistent.
    assert(point.size() == mNumberDimension);

    // Get spatial index.
    kdres *set = kd_nearest(mKdTree, point.data());
    auto spatial_index = *(int *) kd_res_item_data(set);
    kd_res_free(set);

    // Get parameter index.
    int i = 0;
    int parameter_index;
    for (auto &name: mNodalVariableNames) { if (name == parameter_name) parameter_index = i; i++; }

    return mNodalVariables[parameter_index * mNumberNodalVariables + spatial_index];

}
