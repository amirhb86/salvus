#include <salvus.h>
#include <petsc.h>

#include "hdf5.h"
#include "hdf5_hl.h"

void Options::setOptions() {

  /* Error helpers. */
  std::string epre = "Critical option ";
  std::string epst = " not set. Exiting.";

  /* To hold options. */
  PetscInt  int_buffer;
  PetscBool parameter_set, testing, static_problem;
  PetscReal real_buffer;
  char      char_buffer[PETSC_MAX_PATH_LEN];


  /* Set this so that options don't fail if we're unit testing. */
  PetscOptionsGetBool(NULL, NULL, "--testing", &testing, &parameter_set);
  if (parameter_set) {
    testing = PETSC_TRUE;
  } else {
    testing = PETSC_FALSE;
  }
  /* Set this so that options don't fail if we consider a static (i.e., not time-dependent) problem. */
  PetscOptionsGetBool(NULL, NULL, "--static-problem", &static_problem, &parameter_set);
  if (parameter_set) {
    static_problem = PETSC_TRUE;
  } else {
    static_problem = PETSC_FALSE;
  }

  /********************************************************************************
                          Spacial discretization and model.
  ********************************************************************************/
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
  /* TODO: Get this from the problem itself. */
  PetscOptionsGetInt(NULL, NULL, "--dimension", &int_buffer, &parameter_set);
  if (parameter_set) {
    mNumDim = int_buffer;
  } else {
    if (! testing) throw std::runtime_error(epre + "--dimension" + epst);
  }

  /********************************************************************************
                              Time-dependent problems.
  ********************************************************************************/
  PetscOptionsGetReal(NULL, NULL, "--duration", &real_buffer, &parameter_set);
  if (parameter_set) {
    mDuration = real_buffer;
  }
  else {
    if (! testing && ! static_problem ) throw std::runtime_error(epre + "--duration" + epst);
  }
  PetscOptionsGetReal(NULL, NULL, "--time-step", &real_buffer, &parameter_set);
  if (parameter_set) {
    mTimeStep = real_buffer;
  } else {
    if (! testing && ! static_problem ) throw std::runtime_error(epre + "--time-step" + epst);
  }
  // compute number of time steps and ensure that it is integer
  // (i.e., adjust (decrease) mTimeStep if necessary)
  if (! testing && ! static_problem ) {
    mNumTimeSteps = std::ceil(mDuration / mTimeStep);
    mTimeStep = mDuration / (double)(mNumTimeSteps);
    std::cout << "mDuration " << mDuration << ", mTimeStep " << mTimeStep << ", mNumTimeSteps " << mNumTimeSteps << std::endl;
  }
  else {
    mNumTimeSteps = 0;
  }


  /********************************************************************************
                                     Boundaries.
  ********************************************************************************/
  char *bounds[PETSC_MAX_PATH_LEN]; PetscInt num_bnd = PETSC_MAX_PATH_LEN;
  PetscOptionsGetStringArray(NULL, NULL, "--homogeneous-dirichlet", bounds, &num_bnd, &parameter_set);
  if (parameter_set) {
    for (PetscInt i = 0; i < num_bnd; i++) { mHomogeneousDirichletBoundaries.push_back(bounds[i]); }
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
    PetscOptionsGetString(NULL, NULL, "--movie-file-name", char_buffer, PETSC_MAX_PATH_LEN, &parameter_set);
    if (parameter_set) {
      mMovieFile = std::string(char_buffer);
    } else {
      if (! testing)
        throw std::runtime_error("Movie requested, but no output file specified."
                                     " Set --movie-file-name. Exiting.");
    }
    PetscOptionsGetString(NULL, NULL, "--movie-field", char_buffer, PETSC_MAX_PATH_LEN, &parameter_set);
    if (parameter_set) {
      mMovieFields.emplace_back(char_buffer);
    } else {
      if (! testing)
        throw std::runtime_error("Movie requested, but no fields were specified."
                                     " Set --movie-field. Exiting.");
    }
    PetscOptionsGetInt(NULL, NULL, "--save-frame-every", &mSaveFrameEvery, &parameter_set);
    if (!parameter_set) { mSaveFrameEvery = 1; }
  }

  /* TODO: Add options so that sources can be parsed from a file. */
  /********************************************************************************
                                    Sources.
  ********************************************************************************/
  PetscOptionsGetString(NULL, NULL, "--source-file-name", char_buffer, PETSC_MAX_PATH_LEN, &parameter_set);
  if (parameter_set) {
    mSourceFileName = std::string(char_buffer);
    hid_t           file;
    herr_t          status;
    H5G_info_t      info_buffer;

    file = H5Fopen (mSourceFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    status = H5Gget_info (file, &info_buffer);

    mNumSrc = info_buffer.nlinks;

    H5LTget_attribute_string (file, "/", "type", char_buffer);
    mSourceType = std::string(char_buffer);
        
    for ( auto i = 0; i < info_buffer.nlinks; i++) {
      H5Gget_objname_by_idx(file, i, char_buffer, PETSC_MAX_PATH_LEN ) ;
      mSourceNames.push_back(std::string(char_buffer));
      
      std::vector<double> loc(mNumDim);
      H5LTget_attribute_double(file, mSourceNames.at(i).c_str(), "location", loc.data());
      mSrcLocX.push_back(loc[0]);
      mSrcLocY.push_back(loc[1]);
      if (mNumDim == 3) { mSrcLocZ.push_back(loc[2]); }

      H5LTget_attribute_int(file, mSourceNames.at(i).c_str(), "num-components", &int_buffer);
      mSrcNumComponents.push_back(int_buffer);

      if (mSourceType == "ricker") {
        double double_buffer;
        int int_buffer;
        H5LTget_attribute_double(file, mSourceNames.at(i).c_str(), "ricker-amplitude", &double_buffer);
        mSrcRickerAmplitude.push_back(double_buffer);
        H5LTget_attribute_double(file, mSourceNames.at(i).c_str(), "ricker-center-freq", &double_buffer);
        mSrcRickerCenterFreq.push_back(double_buffer);
        H5LTget_attribute_double(file, mSourceNames.at(i).c_str(), "ricker-time-delay", &double_buffer);
        mSrcRickerTimeDelay.push_back(double_buffer);
      }
    }
    status = H5Fclose (file);

  } else {
    PetscOptionsGetInt(NULL, NULL, "--number-of-sources", &int_buffer, &parameter_set);
    if (parameter_set) {
      mNumSrc = int_buffer;
    } else {
      mNumSrc = 0;
    }

    if (mNumSrc > 0) {

      mSrcLocX.resize(mNumSrc); mSrcLocY.resize(mNumSrc);

      PetscOptionsGetString(NULL, NULL, "--source-type", char_buffer, PETSC_MAX_PATH_LEN, &parameter_set);
      if (parameter_set) {
        mSourceType = std::string(char_buffer);
      } else {
        if (! testing)
        throw std::runtime_error("Sources were requested, but source type was not specified. "
                                     "Possibilities are: --source-type [ ricker, file ].");
      }

      PetscInt n_par = mNumSrc; std::string err = "Incorrect number of source parameters: ";
      PetscOptionsGetScalarArray(NULL, NULL, "--source-location-x", mSrcLocX.data(), &n_par, NULL);
      if (n_par != mNumSrc) { throw std::runtime_error(err + "x locations"); }
      PetscOptionsGetScalarArray(NULL, NULL, "--source-location-y", mSrcLocY.data(), &n_par, NULL);
      if (n_par != mNumSrc) { throw std::runtime_error(err + "y locations"); }
      if (mNumDim == 3) {
        mSrcLocZ.resize(mNumSrc);
        PetscOptionsGetScalarArray(NULL, NULL, "--source-location-z", mSrcLocZ.data(), &n_par, NULL);
        if (n_par != mNumSrc) { throw std::runtime_error(err + "z locations"); }
      }

      mSrcNumComponents.resize(mNumSrc);
      PetscOptionsGetIntArray(NULL, NULL, "--source-num-components", mSrcNumComponents.data(), &n_par, NULL);
      if (n_par != mNumSrc) { throw std::runtime_error(err + "--source-num-components"); }

      if (mSourceType == "ricker") {
        n_par = mNumSrc;
        mSrcRickerTimeDelay.resize(mNumSrc);
        mSrcRickerAmplitude.resize(mNumSrc);
        mSrcRickerCenterFreq.resize(mNumSrc);  
        PetscOptionsGetScalarArray(NULL, NULL, "--ricker-amplitude", mSrcRickerAmplitude.data(), &n_par, NULL);
        if (n_par != mNumSrc) { throw std::runtime_error(err + "--ricker-amplitude."); }
        PetscOptionsGetScalarArray(NULL, NULL, "--ricker-time-delay", mSrcRickerTimeDelay.data(), &n_par, NULL);
        if (n_par != mNumSrc) { throw std::runtime_error(err + "--ricker-time-delay"); }
        PetscOptionsGetScalarArray(NULL, NULL, "--ricker-center-freq", mSrcRickerCenterFreq.data(), &n_par, NULL);
        if (n_par != mNumSrc) { throw std::runtime_error(err + "--ricker-center-freq"); }
      } else {
        if (! testing)
        throw std::runtime_error("Source type " + mSourceType + " not recognized.");
      }
    }
  }

  /********************************************************************************
                                    Receivers.
  ********************************************************************************/
  PetscOptionsGetInt(NULL, NULL, "--number-of-receivers", &int_buffer, &parameter_set);
  if (parameter_set) {
    mNumRec = int_buffer;
  } else {
    mNumRec = 0;
  }

  if (mNumRec > 0) {

    mRecLocX.resize(mNumRec); mRecLocY.resize(mNumRec);

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
    if (mNumDim == 3) {
      mRecLocZ.resize(mNumRec);
      PetscOptionsGetScalarArray(NULL, NULL, "--receiver-location-z", mRecLocZ.data(), &n_par, NULL);
      if (n_par != mNumRec) { throw std::runtime_error(err + "--receiver-location-z"); }
    }
  }
}
