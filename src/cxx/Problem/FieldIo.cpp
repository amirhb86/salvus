//
// Created by Michael Afanasiev on 2016-03-16.
//

#include <iostream>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "FieldIo.h"
#include "../Utilities/Utilities.h"

FieldIo *FieldIo::factory(Options options) {
    return new FieldIoHdf5(options);
}

FieldIo::FieldIo(Options options) {

    mPolynomialDegree = options.PolynomialOrder();
    mPhysicsSystem = options.PhysicsSystem();

}


FieldIoHdf5::FieldIoHdf5(Options options): FieldIo(options) {
    std::string fname = "/users/michaelafanasiev/Desktop/test.h5";
    mPlistId = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(mPlistId, PETSC_COMM_WORLD, MPI::INFO_NULL);
    mFileId = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, mPlistId);
    H5Pclose(mPlistId);
}

FieldIoHdf5::~FieldIoHdf5() {
    H5Fclose(mFileId);
}

void FieldIoHdf5::save(const Eigen::MatrixXd &val,
                       const double &time) {

    for (int i = 0; i < mPolynomialDegree; i ++) {
        mSnapshot.nByte.push_back(MPI::COMM_WORLD.Get_rank());
        mSnapshot.byte1.push_back(MPI::COMM_WORLD.Get_rank());
        mSnapshot.byte2.push_back(MPI::COMM_WORLD.Get_rank());
        mSnapshot.byte4.push_back(MPI::COMM_WORLD.Get_rank());
        mSnapshot.byte8.push_back(MPI::COMM_WORLD.Get_rank());
        mSnapshot.scale.push_back(MPI::COMM_WORLD.Get_rank());
        mSnapshot.offset.push_back(MPI::COMM_WORLD.Get_rank());
    }
}

void FieldIoHdf5::read(const double &time) {

}

void FieldIoHdf5::_write_disk() {

    // Name this timedump.
    std::string dump_name = "/TIMEDUMP";

    // Num of datasets to write.
    size_t n_levels = 7;

    // Total size of each dataset (over all processors).
    std::vector<hsize_t> total(n_levels, 0);

    // offset value from beginning of parallel file (processor-specific).
    std::vector<hsize_t> offset(n_levels, 0);

    // Size of each individual array (processor-specific).
    std::vector<hsize_t> size {mSnapshot.nByte.size(), mSnapshot.byte1.size(),
                               mSnapshot.byte2.size(), mSnapshot.byte4.size(),
                               mSnapshot.byte8.size(), mSnapshot.scale.size(),
                               mSnapshot.offset.size()};

    // Datatype of each array.
    std::vector<int> type {H5T_NATIVE_CHAR, H5T_NATIVE_CHAR, H5T_NATIVE_SHORT,
                           H5T_NATIVE_FLOAT, H5T_NATIVE_DOUBLE, H5T_NATIVE_FLOAT,
                           H5T_NATIVE_FLOAT};

    // Name of each dataset.
    std::vector<std::string> names {"NBYTE", "BYTE1", "BYTE2", "BYTE4", "BYTE8",
                                    "SCALE", "OFFSET"};

    // Send offsets to all processors.
    MPI::COMM_WORLD.Scan(size.data(), offset.data(), n_levels,
                         MPI::UNSIGNED_LONG_LONG, MPI::SUM);

    // Send total size to all processors.
    MPI::COMM_WORLD.Allreduce(size.data(), total.data(), size.size(),
                              MPI::UNSIGNED_LONG_LONG, MPI::SUM);

    // Write timedump group.
    hid_t group_id = H5Gcreate(mFileId, dump_name.c_str(), H5P_DEFAULT,
                               H5P_DEFAULT, H5P_DEFAULT);

    // Write time attributes to group.
    hsize_t time_size = mSnapshot.time.size();
    hid_t filespace = H5Screate_simple(1, &time_size, NULL);
    hid_t attr_id = H5Acreate(group_id, "TIMESTEPS", H5T_NATIVE_FLOAT, filespace,
                              H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_FLOAT, mSnapshot.time.data());
    H5Sclose(filespace);
    H5Aclose(attr_id);
    H5Gclose(group_id);

    // TODO: IS THERE A BETTER WAY TO DO THIS?
    // Allocate array to save processor-offsets
    std::vector<unsigned long long> write_offsets (MPI::COMM_WORLD.Get_size());

    // Loop through all datasets and write.
    for (int i = 0; i < n_levels; i++) {

        // Create dataspaces.
        std::string fullname = dump_name + "/" + names[i];
        hid_t filespace = H5Screate_simple(1, &total[i], NULL);
        hid_t dset_id = H5Dcreate(mFileId, fullname.c_str(), type[i], filespace,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        // Define memory space on each processor.
        hsize_t num = size[i];
        hsize_t off = offset[i] - num;
        hid_t memspace = H5Screate_simple(1, &num, NULL);

        // Select hyperslab.
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &off, NULL, &num, NULL);

        // Create property list and write file.
        hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // Get processor offsets onto rank 0.
        MPI_Allgather(&off, 1, MPI::UNSIGNED_LONG_LONG, write_offsets.data(),
                      1, MPI::UNSIGNED_LONG_LONG, PETSC_COMM_WORLD);

        // Write processor offsets.
        hsize_t comm_size = MPI::COMM_WORLD.Get_size();
        hid_t attrspace = H5Screate_simple(1, &comm_size, NULL);
        hid_t attr_id = H5Acreate(dset_id, "PROCCESSOR_OFFSETS", H5T_NATIVE_LLONG, attrspace,
                                  H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_id, H5T_NATIVE_LLONG, write_offsets.data());
        H5Sclose(attrspace);
        H5Aclose(attr_id);

        // Write to different types depending on our array.
        switch (i) {

            case 0:
                H5Dwrite(dset_id, type[i], memspace,
                         filespace, plist_id, mSnapshot.nByte.data());
                break;

            case 1:
                H5Dwrite(dset_id, type[i], memspace,
                         filespace, plist_id, mSnapshot.byte1.data());
                break;

            case 2:
                H5Dwrite(dset_id, type[i], memspace,
                         filespace, plist_id, mSnapshot.byte2.data());
                break;

            case 3:
                H5Dwrite(dset_id, type[i], memspace,
                         filespace, plist_id, mSnapshot.byte4.data());
                break;

            case 4:
                H5Dwrite(dset_id, type[i], memspace,
                         filespace, plist_id, mSnapshot.byte8.data());
                break;

            case 5:
                H5Dwrite(dset_id, type[i], memspace,
                         filespace, plist_id, mSnapshot.scale.data());
                break;

            case 6:
                H5Dwrite(dset_id, type[i], memspace,
                         filespace, plist_id, mSnapshot.offset.data());
                break;

            default:
                break;

        }

        // Release resources.
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Pclose(plist_id);
    }

}

void FieldIoHdf5::write(const double &time) {

    mSnapshot.time.push_back(float (time));

    double eps = 1e-6;
    if (time > 0.01-eps) { _write_disk(); }

}
