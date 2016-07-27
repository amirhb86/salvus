#include <salvus.h>
#include <petsc.h>

void Options::setOptions() {

  /* Error helpers. */
  std::string epre = "Critical option ";
  std::string epst = " not set. Exiting.";

  /* To hold options. */
  PetscInt  int_buffer;
  PetscBool parameter_set, testing;
  PetscReal real_buffer;
  char      char_buffer[PETSC_MAX_PATH_LEN];

  try {

    /* Set this so that options don't fail if we're unit testing. */
    PetscOptionsGetBool(NULL, NULL, "--testing", &testing, &parameter_set);
    if (parameter_set) {
      testing = PETSC_TRUE;
    } else {
      testing = PETSC_FALSE;
    }

    /********************************************************************************
                                 Required options.
    ********************************************************************************/
    PetscOptionsGetReal(NULL, NULL, "--duration", &real_buffer, &parameter_set);
    if (parameter_set) {
      mDuration = real_buffer;
    }
    else {
      if (! testing) throw std::runtime_error(epre + "--duration" + epst);
    }
    PetscOptionsGetReal(NULL, NULL, "--time-step", &real_buffer, &parameter_set);
    if (parameter_set) {
      mTimeStep = real_buffer;
    } else {
      if (! testing) throw std::runtime_error(epre + "--time-step" + epst);
    }
    PetscOptionsGetString(NULL, NULL, "--mesh-file", char_buffer, PETSC_MAX_PATH_LEN, &parameter_set);
    if (parameter_set) {
      mMeshFile = std::string(char_buffer);
    } else {
      if (! testing) throw std::runtime_error(epre + "--mesh-file" + epst);
    }
    PetscOptionsGetString(NULL, NULL, "--model-file", char_buffer, PETSC_MAX_PATH_LEN, &parameter_set);
    if (parameter_set) {
      mModelFile = std::string(char_buffer);
    } else {
      if (! testing) throw std::runtime_error(epre + "--model-file" + epst);
    }
    PetscOptionsGetInt(NULL, NULL, "--polynomial-order", &int_buffer, &parameter_set);
    if (parameter_set) {
      mPolynomialOrder = int_buffer;
    } else {
      if (! testing) throw std::runtime_error(epre + "--polynomial-order" + epst);
    }


    /********************************************************************************
                                         Movies.
    ********************************************************************************/
    PetscOptionsGetBool(NULL, NULL, "--save-movie", &mSaveMovie, &parameter_set);
    if (parameter_set) {
      mSaveMovie = PETSC_TRUE;
    } else {
      mSaveMovie = PETSC_FALSE;
    }

    if (mSaveMovie) {
      PetscOptionsGetString(NULL, NULL, "--output-movie-file-name", char_buffer, PETSC_MAX_PATH_LEN, &parameter_set);
      if (parameter_set) {
        mOutputMovieFile = std::string(char_buffer);
      } else {
        if (! testing)
        throw std::runtime_error("Movie requested, but no output file specified."
                                     " Set --output-movie-file-name. Exiting.");
      }
      PetscOptionsGetInt(NULL, NULL, "--save-frame-every", &mSaveFrameEvery, &parameter_set);
      if (!parameter_set) { mSaveFrameEvery = 1; }
    }

    /* TODO: Add options so that sources can be parsed from a file. */
    /********************************************************************************
                                      Sources.
    ********************************************************************************/
    PetscOptionsGetInt(NULL, NULL, "--number-of-sources", &int_buffer, &parameter_set);
    if (parameter_set) {
      mNumSrc = int_buffer;
    } else {
      mNumSrc = 0;
    }

    if (mNumSrc > 0) {

      mSrcLocX.resize(mNumSrc); mSrcLocY.resize(mNumSrc); mSrcLocZ.resize(mNumSrc);

      PetscOptionsGetString(NULL, NULL, "--source-type", char_buffer, PETSC_MAX_PATH_LEN, &parameter_set);
      if (parameter_set) {
        mSourceType = std::string(char_buffer);
      } else {
        throw std::runtime_error("Sources were requested, but source type was not specified. Exiting.");
      }

      PetscInt n_par = mNumSrc; std::string err = "Incorrect number of source parameters: ";
      PetscOptionsGetScalarArray(NULL, NULL, "--source-location-x", mSrcLocX.data(), &n_par, NULL);
      if (n_par != mNumSrc) { throw std::runtime_error(err + "x locations"); }
      PetscOptionsGetScalarArray(NULL, NULL, "--source-location-y", mSrcLocY.data(), &n_par, NULL);
      if (n_par != mNumSrc) { throw std::runtime_error(err + "y locations"); }
      PetscOptionsGetScalarArray(NULL, NULL, "--source-location-z", mSrcLocZ.data(), &n_par, NULL);
      if (n_par != mNumSrc && mNumDim == 3) { throw std::runtime_error(err + "z locations"); }

      if (mSourceType == "ricker") {
        mSrcRickerTimeDelay.resize(mNumSrc);
        mSrcRickerAmplitude.resize(mNumSrc);
        mSrcRickerCenterFreq.resize(mNumSrc);
        PetscOptionsGetScalarArray(NULL, NULL, "--ricker-amplitude", mSrcRickerAmplitude.data(), &n_par, NULL);
        if (n_par != mNumSrc) { throw std::runtime_error(err + "ricker-amplitude."); }
        PetscOptionsGetScalarArray(NULL, NULL, "--ricker-time-delay", mSrcRickerTimeDelay.data(), &n_par, NULL);
        if (n_par != mNumSrc) { throw std::runtime_error(err + "ricker-time-delay"); }
        PetscOptionsGetScalarArray(NULL, NULL, "--ricker-center-freq", mSrcRickerCenterFreq.data(), &n_par, NULL);
        if (n_par != mNumSrc) { throw std::runtime_error(err + "ricker-center-freq"); }
      } else {
        throw std::runtime_error("Source type " + mSourceType + " not recognized.");
      }
    }

    /********************************************************************************
                                      Recievers.
    ********************************************************************************/
    PetscOptionsGetInt(NULL, NULL, "--number-of-receivers", &int_buffer, &parameter_set);
    if (parameter_set) {
      mNumRec = int_buffer;
    } else {
      mNumRec = 0;
    }

    if (mNumRec > 0) {

      mRecLocX.resize(mNumRec); mRecLocY.resize(mNumRec); mRecLocZ.resize(mNumRec);

      PetscOptionsGetString(NULL, NULL, "--receiver-file-name", char_buffer, PETSC_MAX_PATH_LEN, &parameter_set);
      if (parameter_set) {
        mReceiverFileName = std::string(char_buffer);
      } else {
        if (! testing)
        throw std::runtime_error("Receivers were requested, but no output file was specfied.");
      }

      char *names[PETSC_MAX_PATH_LEN];
      PetscInt n_par = mNumRec; std::string err = "Incorrect number of reciever parameters: ";
      PetscOptionsGetStringArray(NULL, NULL, "--receiver-names", names, &n_par, NULL);
      for (PetscInt i = 0; i < n_par; i++) { mRecNames.push_back(names[i]); }
      if (n_par != mNumRec) { throw std::runtime_error(err + "--receiver-names"); }
      PetscOptionsGetScalarArray(NULL, NULL, "--receiver-location-x", mRecLocX.data(), &n_par, NULL);
      if (n_par != mNumRec) { throw std::runtime_error(err + "--receiver-location-x"); }
      PetscOptionsGetScalarArray(NULL, NULL, "--receiver-location-y", mRecLocY.data(), &n_par, NULL);
      if (n_par != mNumRec) { throw std::runtime_error(err + "--receiver-location-y"); }
      PetscOptionsGetScalarArray(NULL, NULL, "--receiver-location-z", mRecLocZ.data(), &n_par, NULL);
      if (n_par != mNumRec && mNumDim == 3) { throw std::runtime_error(err + "--receiver-location-z"); }
    }

  } catch (std::exception &e) {

    /* It's ok to screw up the options during a test. Otherwise, propagate error. */
    if (testing) {
      throw e;
    } else {
      LOG() << e.what();
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
  }

}
