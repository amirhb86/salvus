#include <petscsys.h>
#include <Utilities/Options.h>
#include <Utilities/Utilities.h>
#include <Utilities/Logging.h>

PetscErrorCode Options::setOptions() {

  // TODO!! DIMENSION ONLY 2 FOR NOW!!
  mDimension = 2;

  PetscInt int_buffer;
  PetscBool parameter_set;
  double real_buffer;

  // Set this so that options don't fail if we're unit testing.
  PetscBool testing;
  PetscOptionsGetBool(NULL, "--testing", &testing, &parameter_set);
  if (! parameter_set) { testing = PETSC_FALSE; }

  PetscOptionsGetReal(NULL, "--duration", &real_buffer, &parameter_set);
  if (parameter_set) { mDuration = real_buffer; }
  else {
    if (! testing)
    { SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Duration must be given via --duration"); }
  }

  PetscOptionsGetReal(NULL, "--time_step", &real_buffer, &parameter_set);
  if (parameter_set) { mTimeStep = real_buffer; }
  else {
    // time step will be set automatically
    mTimeStep = -1;
  }

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
  if (parameter_set) { mNumberSources = int_buffer; } else { mNumberSources = 0; }

  if (mNumberSources > 0) {
    mSourceLocationX.resize(mNumberSources);
    mSourceLocationY.resize(mNumberSources);
    mSourceLocationZ.resize(mNumberSources);
    mSourceRickerTimeDelay.resize(mNumberSources);
    mSourceRickerAmplitude.resize(mNumberSources);
    mSourceRickerCenterFreq.resize(mNumberSources);

    PetscOptionsGetString(NULL, "--source_type", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mSourceType = std::string(char_buffer); }

    PetscOptionsGetScalarArray(NULL, "--source_location_x", mSourceLocationX.data(), &mNumberSources, NULL);
    PetscOptionsGetScalarArray(NULL, "--source_location_y", mSourceLocationY.data(), &mNumberSources, NULL);
    PetscOptionsGetScalarArray(NULL, "--source_location_z", mSourceLocationZ.data(), &mNumberSources, NULL);
    PetscOptionsGetScalarArray(NULL, "--ricker_amplitude", mSourceRickerAmplitude.data(), &mNumberSources, NULL);
    PetscOptionsGetScalarArray(NULL, "--ricker_time_delay", mSourceRickerTimeDelay.data(), &mNumberSources, NULL);
    PetscOptionsGetScalarArray(NULL, "--ricker_center_freq", mSourceRickerCenterFreq.data(), &mNumberSources, NULL);
  }

  // Receivers.
  PetscOptionsGetInt(NULL, "--number_of_receivers", &int_buffer, &parameter_set);
  if (parameter_set) { mNumRec = int_buffer; } else { mNumRec = 0; }
  if (mNumRec > 0) {
    PetscOptionsGetString(NULL, "--receiver_file_name", char_buffer, PETSC_MAX_PATH_LEN,
                          &parameter_set);
    if (parameter_set) { mReceiverFileName = std::string(char_buffer); }
    mRecLocX1.resize(mNumRec);
    mRecLocX2.resize(mNumRec);
    if (mDimension == 3) {
      mRecLocX3.resize(mNumRec);
    }

    int num_max_recs = 1000;
    char *maxrecs[num_max_recs];
    PetscOptionsGetStringArray(NULL, "--receiver_names", maxrecs, &num_max_recs, NULL);
    PetscOptionsGetScalarArray(NULL, "--receiver_location_x1", mRecLocX1.data(), &mNumRec, NULL);
    PetscOptionsGetScalarArray(NULL, "--receiver_location_x2", mRecLocX2.data(), &mNumRec, NULL);
    if (mDimension == 3) {
      PetscOptionsGetScalarArray(NULL, "--receiver_location_x3", mRecLocX3.data(), &mNumRec, NULL);
    }

    for (int i = 0; i < num_max_recs; i++) { mRecNames.push_back(maxrecs[i]); }
  }

  int num_dirichlet_boundaries = 256;
  char *dirichlet_boundaries[256];
  PetscOptionsGetStringArray(NULL,
                             "--dirichlet-boundaries",
                             dirichlet_boundaries,
                             &num_dirichlet_boundaries,
                             &parameter_set);
  if (parameter_set) {
    printf("Using following for dirichlet boundaries:");
    for (int i = 0; i < num_dirichlet_boundaries; i++) {
      printf("%s,", dirichlet_boundaries[i]);
      mDirichletBoundaryNames.push_back(dirichlet_boundaries[i]);
    }
    printf("\n");
  }
  else {
    // default
    mDirichletBoundaryNames.push_back("dirichlet");
  }

  // parameters for movies (to save or not, and how often)
  PetscOptionsGetBool(NULL, "--saveMovie", &mSaveMovie, &parameter_set);
  if (!parameter_set) { mSaveMovie = PETSC_FALSE; }
  PetscOptionsGetInt(NULL, "--saveFrameEvery", &mSaveFrameEvery, &parameter_set);
  if (!parameter_set) { mSaveFrameEvery = 1; }

  // parameters for movies (to save or not, and how often)
  PetscOptionsGetBool(NULL, "--displayDiagnostics", &mDisplayDiagnostics, &parameter_set);
  if (!parameter_set) { mDisplayDiagnostics = PETSC_TRUE; }
  PetscOptionsGetInt(NULL, "--displayDiagnosticsEvery", &mDisplayDiagnosticsEvery, &parameter_set);
  if (!parameter_set) { mDisplayDiagnosticsEvery = 10; }
  
  // for testing ICs against exact solution
  PetscOptionsGetBool(NULL, "--testIC", &mTestIC, &parameter_set);
  if (!parameter_set) { mTestIC = PETSC_FALSE; }
  if (mTestIC) {
    // these will issue unused parameter warning if testIC is false
    PetscOptionsGetReal(NULL, "--IC-center-x", &real_buffer, &parameter_set);
    if (parameter_set) { mCenter_x = real_buffer; }
    PetscOptionsGetReal(NULL, "--IC-center-y", &real_buffer, &parameter_set);
    if (parameter_set) { mCenter_y = real_buffer; }
    PetscOptionsGetReal(NULL, "--IC-center-z", &real_buffer, &parameter_set);
    if (parameter_set) { mCenter_z = real_buffer; }
    PetscOptionsGetReal(NULL, "--IC-square-side-L", &real_buffer, &parameter_set);
    if (parameter_set) { mSquareSide_L = real_buffer; }
  }

  // MAKE THESE COMMAND LINE OPTIONS EVENTUALLY.
  mTimeStepType = "newmark";

  if (mMeshType == "newmark") {
    if (mTestIC) mProblemType = "newmark_testing";
    else mProblemType = "newmark_general";
  }

  // No error
  return 0;
}

void Options::SetTimeStep(double timestep) {
  
    double global_timestep;
    MPI_Allreduce(&timestep,&global_timestep,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    DEBUG() << "Local timestep " << timestep << " vs. global timestep " << global_timestep;
    mTimeStep = global_timestep;
    
}
