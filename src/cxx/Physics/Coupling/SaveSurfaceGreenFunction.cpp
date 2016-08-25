#include <mpi.h>
#include <Utilities/Options.h>
#include <Utilities/Utilities.h>
#include <Physics/SaveSurfaceGreenFunction.h>

template <typename Element>
MPI_Group SaveSurfaceGreenFunction<Element>::mMpiGroup = NULL;
template <typename Element>
MPI_Comm SaveSurfaceGreenFunction<Element>::mMpiComm = NULL;
template <typename Element>
hid_t SaveSurfaceGreenFunction<Element>::mFileId = 0;
template <typename Element>
PetscInt SaveSurfaceGreenFunction<Element>::mNumTimeStep;
template <typename Element>
PetscInt SaveSurfaceGreenFunction<Element>::mNumMixins = 0;
template <typename Element>
UncompressedWavefieldContainer<PetscReal> SaveSurfaceGreenFunction<Element>::container;

template <typename Element>
SaveSurfaceGreenFunction<Element>::SaveSurfaceGreenFunction(
    std::unique_ptr<Options> const &options) : Element(options) {

  if (options->NumTimeSteps() <= 0) {
    throw std::runtime_error("Invalid number of timesteps provided to SaveSurfaceGreenFunction. "
                                 "\nNumber of timesteps: " + std::to_string(mNumTimeStep));
  }

  mNumTimeStep = options->NumTimeSteps();
  mCountTimeSteps = 0;
  mNumComponents = Element::PullElementalFields().size();
  mMixNum = mNumMixins++;

}

template <typename Element>
SaveSurfaceGreenFunction<Element>::~SaveSurfaceGreenFunction() {
  write();
  mNumMixins--;
}

template <typename Element>
void SaveSurfaceGreenFunction<Element>::write() {

  if (!mFileId) {

    if (!container.size()) return;

    /* Create communicator from ranks. */
    std::vector<PetscInt> ranks = utilities::GetWorldRanksForTag("SaveSurface");
    MPI_Group world_group;
    MPI_Comm_group(PETSC_COMM_WORLD, &world_group);
    MPI_Group_incl(world_group, ranks.size(), ranks.data(), &mMpiGroup);
    MPI_Comm_create_group(PETSC_COMM_WORLD, mMpiGroup, 0, &mMpiComm);

    /* Get total number of elements. */
    PetscInt num_all_mixins;
    MPI_Allreduce(&mNumMixins, &num_all_mixins, 1, MPIU_INT, MPI_SUM, mMpiComm);
    int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    /* Get element offsets. */
    PetscInt mixin_offest;
    MPI_Scan(&mNumMixins, &mixin_offest, 1, MPIU_INT, MPI_SUM, mMpiComm);

    /* Create file with communicator. */
    std::string fname = "./test.h5";
    {
      hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(plist_id, mMpiComm, MPI_INFO_NULL);
      mFileId = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
      H5Pclose(plist_id);
    }

    const hsize_t dim = 4;
    const hsize_t size[dim] = {(hsize_t) container.tsp(), (hsize_t) num_all_mixins,
                              (hsize_t) container.cmp(), (hsize_t) container.pnt()};

    {
      hid_t filespace = H5Screate_simple(dim, size, NULL);
      hid_t dset_id = H5Dcreate(mFileId, "wavefield_data", container.Hdf5Datatype(),
                                filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      H5DSset_scale(dset_id, Element::Name().c_str());
      H5DSset_label(dset_id, 0, "time_step");
      H5DSset_label(dset_id, 1, "element");
      std::string components = "[";
      for (auto &f: Element::PullElementalFields()) { components += " " + f + " "; }
      components += "]";
      H5DSset_label(dset_id, 2, components.c_str());
      H5DSset_label(dset_id, 3, "point");

      H5Sclose(filespace);
      const hsize_t off[dim] = {0, (hsize_t) mixin_offest-mNumMixins, 0, 0};
      filespace = H5Dget_space(dset_id);
      const hsize_t size_this_proc[dim] = {size[0], (hsize_t) mNumMixins, size[2], size[3] };
      hid_t memspace = H5Screate_simple(dim, size_this_proc, NULL);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, off, NULL, size_this_proc, NULL);
      hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Dwrite(dset_id, container.Hdf5Datatype(), memspace, filespace, plist_id, &container.data());
      H5Dclose(dset_id);
      H5Sclose(filespace);
      H5Sclose(memspace);
      H5Pclose(plist_id);

    }

    H5Fclose(mFileId);
  }

}

template <typename Element>
void SaveSurfaceGreenFunction<Element>::recordDynamicFields
    (const Eigen::Ref<const RealMat>& field) {

  /* Initialize the container only once. */
  if (!container.size() && mMixNum == 0) {
    container.resize(mNumTimeStep, mNumMixins, mNumComponents, Element::NumIntPnt());
  }

  /* Save the component. */
  for (PetscInt i = 0; i < mNumComponents; ++i) {
    for (PetscInt j = 0; j < Element::NumIntPnt(); ++j) {
      container(mCountTimeSteps, mMixNum, i, j) = field(j, i);
    }
  }
  mCountTimeSteps++;

}


#include <Physics/Scalar.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>

template class SaveSurfaceGreenFunction<
    SaveSurfaceGreenFunctionTestStub>;

template class SaveSurfaceGreenFunction<
    Scalar<
        TensorQuad<
            QuadP1>>>;
