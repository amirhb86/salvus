#include "QuadNew.h"

class Physics;

template <typename Physics>
void QuadNew<Physics>::prepareStiffness() {
  std::cout << "Prepare Stiffness." << std::endl;
}

template <typename Physics>
void QuadNew<Physics>::attachVertexCoordinates(DM &distributed_mesh) {

  Vec coordinates_local;
  PetscInt coordinate_buffer_size;
  PetscSection coordinate_section;
  PetscReal *coordinates_buffer = NULL;

  DMGetCoordinatesLocal(distributed_mesh, &coordinates_local);
  DMGetCoordinateSection(distributed_mesh, &coordinate_section);
  DMPlexVecGetClosure(distributed_mesh, coordinate_section, coordinates_local, mElmNum,
                      &coordinate_buffer_size, &coordinates_buffer);
  std::vector<PetscReal> coordinates_element(coordinates_buffer, coordinates_buffer + coordinate_buffer_size);
  DMPlexVecRestoreClosure(distributed_mesh, coordinate_section, coordinates_local, mElmNum,
                          &coordinate_buffer_size, &coordinates_buffer);

  for (int i = 0; i < mNumVtx; i++) {
    mVtxCrd(i,0) = coordinates_element[mNumDim * i + 0];
    mVtxCrd(i,1) = coordinates_element[mNumDim * i + 1];
  }

  // Save element center
  mElmCtr << mVtxCrd.col(0).mean(),
      mVtxCrd.col(1).mean();

}

template class QuadNew<AcousticNew>;
