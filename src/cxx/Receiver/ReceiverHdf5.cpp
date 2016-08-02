#include <Utilities/Options.h>
#include <Utilities/Utilities.h>
#include <Receiver/ReceiverHdf5.h>

// Initialize hdf5 static.
hid_t ReceiverHdf5::mFileId;
std::vector<std::string> ReceiverHdf5::mWriteRegisteredFields;

ReceiverHdf5::ReceiverHdf5(std::unique_ptr<Options> const &options) : Receiver(options) {

  // Only create one hdf5 file for all receivers.
  if ((Num() == 0) && (options->ReceiverFileName().size())) {

    // Create file and set access.
    std::string fname = options->ReceiverFileName();
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, PETSC_COMM_WORLD, MPI_INFO_NULL);
    mFileId = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

  }

}

ReceiverHdf5::~ReceiverHdf5() {
  if ((Num() == 0) && mFileId) {
    H5Fclose(mFileId);
  }
}

void ReceiverHdf5::write() {

  // Send maximum length to all processors. Also send all field names (just once).
  hsize_t loc_size = (store.size()) ? store.begin()->second.size() : 0;
  hsize_t n_samples = 0;
  MPI_Allreduce(&loc_size, &n_samples, 1, MPI_UNSIGNED_LONG_LONG,
                MPI_SUM, MPI_COMM_WORLD);

  /* TODO: POSSIBLE BUG! Will only write fields named on receiver zero. */
  if (Num() == 0) {
    int name_rank = 0;
    std::vector<std::string> name;
    for (auto &dict: store) {
      name.push_back(dict.first);
      MPI_Comm_rank(PETSC_COMM_WORLD, &name_rank);
    }
    MPI_Allreduce(MPI_IN_PLACE, &name_rank, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    mWriteRegisteredFields = utilities::broadcastStringVecFromRank(name, name_rank);
  }

  hid_t gid = H5Gcreate(mFileId, Name().c_str(), H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(gid);

  for (auto &field: mWriteRegisteredFields) {

    std::string fullname = "/" + Name() + "/" + field;
    hid_t filespace = H5Screate_simple(1, &n_samples, NULL);
    hid_t dset_id = H5Dcreate(mFileId, fullname.c_str(), H5T_NATIVE_FLOAT, filespace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    // Define memory space on each processor.
    hsize_t num = loc_size;
    hid_t memspace = H5Screate_simple(1, &num, NULL);

    // Select hyperslab.
    hsize_t off = 0;
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &off, NULL, &num, NULL);
    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write if field is present on this processor. Else write zeros. */
    if (store.find(field) != store.end()) {
      H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, store[field].data());
    } else {
      Eigen::VectorXf zeros = Eigen::VectorXf::Zero(store.begin()->second.size());
      H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, zeros.data());
    }
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
  }

}



