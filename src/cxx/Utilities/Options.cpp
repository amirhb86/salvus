//
// Created by Michael Afanasiev on 2016-01-29.
//

#include <petscsys.h>
#include "Options.h"
#include "Utilities.h"


PetscErrorCode Options::setOptions() {

    PetscInt int_buffer;
    PetscBool parameter_set;
    double real_buffer;

    PetscOptionsGetReal(NULL, "--duration", &real_buffer, &parameter_set);
    if (parameter_set) { mDuration = real_buffer; }
    else { SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Duration must be given via --duration");}

    PetscOptionsGetReal(NULL, "--time_step", &real_buffer, &parameter_set);
    if (parameter_set) { mTimeStep = real_buffer; }
    else { SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Time step must be given via --time_step");}

    // String options.O
    char char_buffer[PETSC_MAX_PATH_LEN];

    PetscOptionsGetString(NULL, "--exodus_file_name", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mExodusMeshFile = std::string(char_buffer); }

    PetscOptionsGetString(NULL, "--exodus_model_file_name", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mExodusModelFile = std::string(char_buffer); }

    PetscOptionsGetString(NULL, "--output_movie_file_name", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mOutputMovieFile = std::string(char_buffer); }

    PetscOptionsGetString(NULL, "--mesh_type", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mMeshType = std::string(char_buffer); }

    PetscOptionsGetString(NULL, "--element_shape", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mElementShape = std::string(char_buffer); }

    PetscOptionsGetString(NULL, "--physics_system", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mPhysicsSystem = std::string(char_buffer); }

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

    if(mNumberSources > 0) {
      PetscOptionsGetString(NULL, "--source_type", char_buffer, PETSC_MAX_PATH_LEN,
                            &parameter_set);
      if (parameter_set) { mSourceType = std::string(char_buffer); }
      
      PetscOptionsGetScalarArray(NULL, "--source_location_x", mSourceLocationX.data(), &mNumberSources, NULL);
      //    PetscOptionsGetScalarArray(NULL, "--source_location_y", mSourceLocationY.data(), &mNumberSources, NULL);
      PetscOptionsGetScalarArray(NULL, "--source_location_z", mSourceLocationZ.data(), &mNumberSources, NULL);
      PetscOptionsGetScalarArray(NULL, "--ricker_amplitude", mSourceRickerAmplitude.data(), &mNumberSources, NULL);
      PetscOptionsGetScalarArray(NULL, "--ricker_time_delay", mSourceRickerTimeDelay.data(), &mNumberSources, NULL);
      PetscOptionsGetScalarArray(NULL, "--ricker_center_freq", mSourceRickerCenterFreq.data(), &mNumberSources, NULL);
    }
    
    // parameters for movies (to save or not, and how often)
    PetscOptionsGetBool(NULL, "--saveMovie", &mSaveMovie,&parameter_set);
    if(!parameter_set) { mSaveMovie = PETSC_FALSE; }
    PetscOptionsGetInt(NULL, "--saveFrameEvery", &mSaveFrameEvery,&parameter_set);
    if(!parameter_set) { mSaveFrameEvery = 1; }
    
    // MAKE THESE COMMAND LINE OPTIONS EVENTUALLY.
    mDimension = 2;
    mTimeStepType = "newmark";

    if(mMeshType == "newmark")
        { mProblemType = "newmark_general"; }

    // No error
    return 0;
}
