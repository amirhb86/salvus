#include <iostream>
#include "ReceiverHdf5.h"

// Initialize hdf5 static.
hid_t ReceiverHdf5::mPlistId;
hid_t ReceiverHdf5::mFileId;
std::vector<std::string> ReceiverHdf5::names;

ReceiverHdf5::ReceiverHdf5(Options options) : Receiver(options) {

  // Only create one hdf5 file for all receivers.
  if (Num() == 0) {

    // Create file and set access.
    std::string fname = "/users/michaelafanasiev/Desktop/test.h5";
    mPlistId = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(mPlistId, PETSC_COMM_WORLD, MPI_INFO_NULL);
    mFileId = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, mPlistId);
    H5Pclose(mPlistId);

  }

  // Save all names globally (for HDF5 dataset creation).
  names.push_back(options.RecNames()[Num()]);

}

ReceiverHdf5::~ReceiverHdf5() {
  if (Num() == 0) {
    H5Fclose(mFileId);
  }
}

void ReceiverHdf5::write() {


}



