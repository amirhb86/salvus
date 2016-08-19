#include <mpi.h>
#include <Source/Source.h>
#include <Source/SourceHdf5.h>
#include <Utilities/Options.h>
#include <stdexcept>
#include <iostream>
#include <Utilities/Logging.h>

#include "hdf5.h"
#include "hdf5_hl.h"


SourceHdf5::SourceHdf5(std::unique_ptr<Options> const &options): Source(options) {

  /* Set locations. */


  SetLocX(options->SrcLocX()[Num()]);
  SetLocY(options->SrcLocY()[Num()]);
  if (options->SrcLocZ().size()) {
    SetLocZ(options->SrcLocZ()[Num()]);
  }

  mSourceFileName = options->SourceFileName();
  mSourceName = options->SrcName()[Num()];
  mNumTimeSteps = options->NumTimeSteps();
  mNumComponents = options->SrcNumComponents()[Num()];
}

Eigen::VectorXd SourceHdf5::fire(const double &time, const PetscInt &time_idx) {

  Eigen::VectorXd src = source_time_function.row(time_idx);
  return src;

}


void SourceHdf5::loadData() {

  hid_t           file;
  herr_t          status;

  file = H5Fopen (mSourceFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  
  std::string dataset_name = mSourceName + "/data";
  source_time_function.setZero(mNumTimeSteps,mNumComponents);
  
  H5LTread_dataset_double(file, dataset_name.c_str(), source_time_function.data());
  
  H5Fclose (file);
}