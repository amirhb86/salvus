//
// Created by Michael Afanasiev on 2016-01-29.
//

#include <petscsys.h>
#include "Options.h"
#include "Utilities.h"


void Options::setOptions() {

    PetscInt int_buffer;
    PetscBool parameter_set;
    char char_buffer[PETSC_MAX_PATH_LEN];

    PetscOptionsGetString(NULL, "--exodus_file_name", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mExodusMeshFile = std::string(char_buffer); }

    PetscOptionsGetString(NULL, "--exodus_model_file_name", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mExodusModelFile = std::string(char_buffer); }

    PetscOptionsGetString(NULL, "--mesh_type", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mMeshType = std::string(char_buffer); }

    PetscOptionsGetString(NULL, "--element_shape", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mElementShape = std::string(char_buffer); }

    PetscOptionsGetString(NULL, "--physics_system", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mPhysicsSystem = std::string(char_buffer); }

    PetscOptionsGetString(NULL, "--source_type", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mSourceType = std::string(char_buffer); }

    PetscOptionsGetInt(NULL, "--polynomial_order", &int_buffer, &parameter_set);
    if (parameter_set) { mPolynomialOrder = int_buffer; }

    // Sources.

    PetscOptionsGetInt(NULL, "--number_of_sources", &int_buffer, &parameter_set);
    if (parameter_set) { mNumberSources = int_buffer; }

    mSourceLocationX.resize(mNumberSources);
    mSourceLocationY.resize(mNumberSources);
    mSourceLocationZ.resize(mNumberSources);
    mSourceRickerTimeDelay.resize(mNumberSources);
    mSourceRickerAmplitude.resize(mNumberSources);
    mSourceRickerCenterFreq.resize(mNumberSources);

    PetscOptionsGetScalarArray(NULL, "--source_location_x", mSourceLocationX.data(), &mNumberSources, NULL);
//    PetscOptionsGetScalarArray(NULL, "--source_location_y", mSourceLocationY.data(), &mNumberSources, NULL);
    PetscOptionsGetScalarArray(NULL, "--source_location_z", mSourceLocationZ.data(), &mNumberSources, NULL);
    PetscOptionsGetScalarArray(NULL, "--ricker_amplitude", mSourceRickerAmplitude.data(), &mNumberSources, NULL);
    PetscOptionsGetScalarArray(NULL, "--ricker_time_delay", mSourceRickerTimeDelay.data(), &mNumberSources, NULL);
    PetscOptionsGetScalarArray(NULL, "--ricker_center_freq", mSourceRickerCenterFreq.data(), &mNumberSources, NULL);

}
