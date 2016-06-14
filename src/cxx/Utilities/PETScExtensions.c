#include <Utilities/PETScExtensions.h>

/* closure helpers */
#undef __FUNCT__
#define __FUNCT__ "_DMPlexVecGetClosure_Depth1_Static_GetIndices"
PetscErrorCode _DMPlexVecGetClosure_Depth1_Static_GetIndices(DM dm, PetscSection section, PetscInt point, PetscInt *csize, PetscInt *idx[])
{
  PetscInt       *ilist;
  const PetscInt *cone, *coneO;
  PetscInt        pStart, pEnd, p, numPoints, size = 0, offset = 0;
  PetscErrorCode  ierr;

  ierr = PetscSectionGetChart(section, &pStart, &pEnd);CHKERRQ(ierr);
  ierr = DMPlexGetConeSize(dm, point, &numPoints);CHKERRQ(ierr);
  ierr = DMPlexGetCone(dm, point, &cone);CHKERRQ(ierr);
  ierr = DMPlexGetConeOrientation(dm, point, &coneO);CHKERRQ(ierr);
  {
    if ((point >= pStart) && (point < pEnd)) {
      PetscInt dof;

      ierr = PetscSectionGetDof(section, point, &dof);CHKERRQ(ierr);
      size += dof;
    }
    for (p = 0; p < numPoints; ++p) {
      const PetscInt cp = cone[p];
      PetscInt       dof;

      if ((cp < pStart) || (cp >= pEnd)) continue;
      ierr = PetscSectionGetDof(section, cp, &dof);CHKERRQ(ierr);
      size += dof;
    }

    PetscMalloc1(size,&ilist);
  }

  size = 0;
  if ((point >= pStart) && (point < pEnd)) {
    PetscInt     dof, off, d;

    ierr = PetscSectionGetDof(section, point, &dof);CHKERRQ(ierr);
    ierr = PetscSectionGetOffset(section, point, &off);CHKERRQ(ierr);
    for (d = 0; d < dof; ++d, ++offset) {
      ilist[offset] = off + d;
    }
    size += dof;
  }
  for (p = 0; p < numPoints; ++p) {
    const PetscInt cp = cone[p];
    PetscInt       o  = coneO[p];
    PetscInt       dof, off, d;
    PetscScalar   *varr;

    if ((cp < pStart) || (cp >= pEnd)) continue;
    ierr = PetscSectionGetDof(section, cp, &dof);CHKERRQ(ierr);
    ierr = PetscSectionGetOffset(section, cp, &off);CHKERRQ(ierr);
    if (o >= 0) {
      for (d = 0; d < dof; ++d, ++offset) {
        ilist[offset] = off + d;
      }
    } else {
      for (d = dof-1; d >= 0; --d, ++offset) {
        ilist[offset] = off + d;
      }
    }
    size += dof;
  }

  *csize = size;
  *idx = ilist;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_DMPlexVecGetClosure_Fields_Static_GetIndices"
PetscErrorCode _DMPlexVecGetClosure_Fields_Static_GetIndices(PetscSection section, PetscInt numPoints, const PetscInt points[], PetscInt numFields, PetscInt *size, PetscInt idx[])
{
  PetscInt       offset = 0, f;
  PetscErrorCode ierr;

  *size = 0;
  for (f = 0; f < numFields; ++f) {
    PetscInt fcomp, p;

    ierr = PetscSectionGetFieldComponents(section, f, &fcomp);CHKERRQ(ierr);

    for (p = 0; p < numPoints*2; p += 2) {
      const PetscInt point = points[p];
      const PetscInt o     = points[p+1];
      PetscInt       fdof, foff, d, c;

      ierr = PetscSectionGetFieldDof(section, point, f, &fdof);CHKERRQ(ierr);
      ierr = PetscSectionGetFieldOffset(section, point, f, &foff);CHKERRQ(ierr);
      if (o >= 0) {
        for (d = 0; d < fdof; ++d, ++offset) idx[offset] = foff + d;
      } else {
        for (d = fdof/fcomp-1; d >= 0; --d) {
          for (c = 0; c < fcomp; ++c, ++offset) {
            idx[offset] = foff + d*fcomp+c;
          }
        }
      }
    }
  }
  *size = offset;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_DMPlexVecGetClosure_Static_GetIndices"
PetscErrorCode _DMPlexVecGetClosure_Static_GetIndices(PetscSection section, PetscInt numPoints, const PetscInt points[], PetscInt *size, PetscInt idx[])
{
  PetscInt       offset = 0, p;
  PetscErrorCode ierr;

  *size = 0;
  for (p = 0; p < numPoints*2; p += 2) {
    const PetscInt point = points[p];
    const PetscInt o     = points[p+1];
    PetscInt       dof, off, d;

    ierr = PetscSectionGetDof(section, point, &dof);CHKERRQ(ierr);
    ierr = PetscSectionGetOffset(section, point, &off);CHKERRQ(ierr);
    if (o >= 0) {
      for (d = 0; d < dof; ++d, ++offset)    idx[offset] = off + d;
    } else {
      for (d = dof-1; d >= 0; --d, ++offset) idx[offset] = off + d;
    }
  }
  *size = offset;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPlexGetClosureIndices"
PetscErrorCode DMPlexGetClosureIndices(DM dm, PetscSection section, PetscInt point, PetscInt *csize, PetscInt *idx[])
{
  PetscSection    clSection;
  IS              clPoints;
  PetscInt       *points = NULL;
  PetscInt       *ilist;
  const PetscInt *clp;
  PetscInt        depth, numFields, numPoints, size;
  PetscErrorCode  ierr;

  if (!csize) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"Must provide a pointer to csize");
  if (!idx) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"Must provide a pointer to idx");
  if (!section) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"Must provide a PetscSection");

  ierr = DMPlexGetDepth(dm, &depth);CHKERRQ(ierr);
  ierr = PetscSectionGetNumFields(section, &numFields);CHKERRQ(ierr);
  if (depth == 1 && numFields < 2) {
    ierr = _DMPlexVecGetClosure_Depth1_Static_GetIndices(dm, section, point, csize, idx);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  /* Get points */
  ierr = PetscSectionGetClosureIndex(section, (PetscObject) dm, &clSection, &clPoints);CHKERRQ(ierr);
  if (!clPoints) {
    PetscInt pStart, pEnd, p, q;

    ierr = PetscSectionGetChart(section, &pStart, &pEnd);CHKERRQ(ierr);
    ierr = DMPlexGetTransitiveClosure(dm, point, PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);
    /* Compress out points not in the section */
    for (p = 0, q = 0; p < numPoints*2; p += 2) {
      if ((points[p] >= pStart) && (points[p] < pEnd)) {
        points[q*2]   = points[p];
        points[q*2+1] = points[p+1];
        ++q;
      }
    }
    numPoints = q;
  } else {
    PetscInt dof, off;

    ierr = PetscSectionGetDof(clSection, point, &dof);CHKERRQ(ierr);
    ierr = PetscSectionGetOffset(clSection, point, &off);CHKERRQ(ierr);
    ierr = ISGetIndices(clPoints, &clp);CHKERRQ(ierr);
    numPoints = dof/2;
    points    = (PetscInt *) &clp[off];
  }

  /* Get array */
  {
    PetscInt asize = 0, dof, p;

    for (p = 0; p < numPoints*2; p += 2) {
      ierr = PetscSectionGetDof(section, points[p], &dof);CHKERRQ(ierr);
      asize += dof;
    }
    *csize = asize;
    PetscMalloc1(asize,&ilist);
  }

  /* Get values  - Reimplemeneted as GetIndices */
  if (numFields > 0) {ierr = _DMPlexVecGetClosure_Fields_Static_GetIndices(section, numPoints, points, numFields, &size, ilist);CHKERRQ(ierr);}
  else               {ierr = _DMPlexVecGetClosure_Static_GetIndices(section, numPoints, points, &size, ilist);CHKERRQ(ierr);}
  /* Cleanup points */
  if (!clPoints) {ierr = DMPlexRestoreTransitiveClosure(dm, point, PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);}
  else           {ierr = ISRestoreIndices(clPoints, &clp);CHKERRQ(ierr);}

  if (size > *csize) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Size of input array %d < actual size %d", *csize, size);

  *csize = size;
  *idx = ilist;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMPlexGetClosureWorkArray"
PetscErrorCode DMPlexGetClosureWorkArray(DM dm,PetscInt npoints,PetscSection section,PetscInt *bsize,PetscScalar **buffer)
{
  PetscErrorCode ierr;
  PetscInt e,max,size;
  PetscScalar *array;
  Vec u;

  if (!section) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"Must provide a valid PetscSection");
  if (!bsize) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"Must provide a pointer to bsize");
  if (!buffer) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"Must provide a pointer to the work buffer");

  ierr = DMGetLocalVector(dm,&u);CHKERRQ(ierr);
  max = PETSC_MIN_INT;
  for (e=0; e<npoints; e++) {
    PetscScalar *vals = NULL;

    ierr = DMPlexVecGetClosure(dm,section,u,e,&size,&vals);CHKERRQ(ierr);
    ierr = DMPlexVecRestoreClosure(dm,section,u,e,&size,&vals);CHKERRQ(ierr);
    if (size > max) { max = size; }
  }
  ierr = DMRestoreLocalVector(dm,&u);CHKERRQ(ierr);

  PetscMalloc1(max,&array);
  *bsize = max;
  *buffer = array;

  PetscFunctionReturn(0);
}


