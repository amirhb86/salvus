#if !defined(__PETSCEXTEN_H)
#define __PETSCEXTEN_H

#include <petsc.h>

PetscErrorCode _DMPlexVecGetClosure_Depth1_Static_GetIndices(
    DM dm, PetscSection section, PetscInt point, PetscInt *csize, PetscInt *idx[]);

PetscErrorCode _DMPlexVecGetClosure_Fields_Static_GetIndices(
    PetscSection section, PetscInt numPoints, const PetscInt points[],
    PetscInt numFields, PetscInt *size, PetscInt idx[]);

PetscErrorCode _DMPlexVecGetClosure_Static_GetIndices(
    PetscSection section, PetscInt numPoints, const PetscInt points[], PetscInt *size,
    PetscInt idx[]);

PetscErrorCode DMPlexGetClosureIndices(DM dm, PetscSection section, PetscInt point,
                                       PetscInt *csize, PetscInt *idx[]);

PetscErrorCode DMPlexGetClosureWorkArray(DM dm, PetscInt npoints,PetscSection section,
                                         PetscInt *bsize, PetscScalar **buffer);

#endif





