#include <petsc.h>
#include <salvus.h>
#include <hdf5.h>
#include "catch.h"

using namespace std;

TEST_CASE("Test receiver functionality", "[receiver]") {

  SECTION("unit") {

    PetscOptionsClear(NULL);
    const char *arg[] =
        {"salvus_test", "--testing", "true",
         "--number-of-receivers", "2",
         "--receiver-file-name", "./test_receiver.h5",
         "--receiver-names", "rec1,rec2",
         "--receiver-location-x", "50000,50000",
         "--receiver-location-y", "50000,90000",
         "--receiver-location-z", "50000,90000",
         NULL};

    char **argv = const_cast<char **> (arg);
    int argc = sizeof(arg) / sizeof(const char *) - 1;
    PetscOptionsInsert(NULL, &argc, &argv, NULL);

    SECTION("general") {

      vector<PetscReal> x{50000, 50000};
      vector<PetscReal> y{50000, 90000};
      vector<PetscReal> z{50000, 90000};

      unique_ptr<Options> options(new Options);
      options->SetDimension(3);
      options->setOptions();

      /* Generate receiver vector. */
      auto receivers = Receiver::Factory(options);

      /* Ensure that the proper number were allocated. */
      REQUIRE(Receiver::NumReceivers() == 2);

      /* Require that all locations were set properly. */
      for (PetscInt i = 0; i < Receiver::NumReceivers(); i++) {
        REQUIRE(receivers[i]->LocX() == x[i]);
        REQUIRE(receivers[i]->LocY() == y[i]);
        REQUIRE(receivers[i]->LocZ() == z[i]);
      }

      /* Test that we're recording something. */
      RealVec trial_vals = RealVec::LinSpaced(10, 1, 10);
      for (PetscInt i = 0; i < trial_vals.size(); i++) {
        receivers[0]->record(trial_vals[i] + 0, "test_0");
        receivers[0]->record(trial_vals[i] + 1, "test_1");
        receivers[1]->record(trial_vals[i] + 2, "test_0");
      }

      /* Write. */
      for (auto &rec: receivers) {
        rec->write();
      }

      /* Require that deletion and deallocation works properly. */
      receivers.back().reset();
      REQUIRE(Receiver::NumReceivers() == 1);
      receivers[0].reset();
      REQUIRE(Receiver::NumReceivers() == 0);

      /* Ensure write was correct. */
      Eigen::VectorXf data(10);
      std::string fname = "test_receiver.h5";
      hid_t fileid = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      hid_t groupid = H5Gopen(fileid, "/rec1", H5P_DEFAULT);
      hid_t dset_id = H5Dopen2(groupid, "/rec1/test_0", H5P_DEFAULT);
      herr_t status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                              data.data());
      H5Dclose(dset_id);
      REQUIRE(data.isApprox(Eigen::VectorXf::LinSpaced(10, 1, 10)));
      dset_id = H5Dopen2(groupid, "/rec1/test_1", H5P_DEFAULT);
      status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                              data.data());
      REQUIRE(data.isApprox(Eigen::VectorXf::LinSpaced(10, 2, 11)));
      H5Dclose(dset_id);
      H5Gclose(groupid);
      groupid = H5Gopen(fileid, "/rec2", H5P_DEFAULT);
      dset_id = H5Dopen2(groupid, "/rec2/test_0", H5P_DEFAULT);
      status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                              data.data());
      REQUIRE(data.isApprox(Eigen::VectorXf::LinSpaced(10, 3, 12)));
      H5Dclose(dset_id);
      dset_id = H5Dopen2(groupid, "/rec2/test_1", H5P_DEFAULT);
      status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data.data());
      REQUIRE(data.isApprox(Eigen::VectorXf::Zero(10)));
      H5Dclose(dset_id);
      H5Gclose(groupid);
      H5Fclose(fileid);

    }

  }

  SECTION("exceptions") {
    unique_ptr<Options> options(new Options);
    options->setOptions();
    options->SetReceiverFileName("test");

    REQUIRE_THROWS_AS(Receiver::Factory(options), std::runtime_error);
  }



}
