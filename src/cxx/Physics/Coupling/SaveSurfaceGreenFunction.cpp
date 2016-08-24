#include <mpi.h>
#include <Utilities/Utilities.h>
#include <Physics/SaveSurfaceGreenFunction.h>

template <typename Element>
MPI_Group SaveSurfaceGreenFunction<Element>::mMpiGroup = NULL;
template <typename Element>
MPI_Comm  SaveSurfaceGreenFunction<Element>::mMpiComm = NULL;

template <typename Element>
SaveSurfaceGreenFunction<Element>::SaveSurfaceGreenFunction(
    std::unique_ptr<Options> const &options) : Element(options) {

}

template <typename Element>
void SaveSurfaceGreenFunction<Element>::recordDynamicFields
    (const Eigen::Ref<const RealMat>& field) {

  if (!mMpiComm) {
    /* Create group from ranks. */
    std::vector<PetscInt> ranks = utilities::GetWorldRanksForTag("SaveSurface");

  }

}

#include <Physics/Scalar.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>

template class SaveSurfaceGreenFunction<
    Scalar<
        TensorQuad<
            QuadP1>>>;
