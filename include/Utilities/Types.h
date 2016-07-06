#pragma once

#include <petsc.h>
#include <Eigen/Dense>
#include <Element/Element.h>

/* Structure managing PETSc vectors. */
struct field {
  field(const std::string name, DM &dm)
  {
    mName = name;
    DMCreateLocalVector(dm, &mLoc);  VecSet(mLoc, 0);
    DMCreateGlobalVector(dm, &mGlb); VecSet(mGlb, 0);
    PetscObjectSetName((PetscObject) mGlb, mName.c_str());
  }
  ~field()          /** < Clean memory */
  {
    if (mGlb) { VecDestroy(&mGlb); }
    if (mLoc) { VecDestroy(&mLoc); }
  }
  std::string mName; /** < Field name (i.e. displacement_x) */
  Vec mGlb;          /** < Global PETSc vector */
  Vec mLoc;          /** < Local PETSc vector */
};

/* Here and some convienience typedefs to save space in functions. */
typedef Eigen::Matrix<PetscInt, Eigen::Dynamic, 1> IntVec;
typedef Eigen::Matrix<PetscReal, Eigen::Dynamic, 1> RealVec;
typedef std::vector<std::unique_ptr<Element>> ElemVec;
typedef std::map<std::string, std::unique_ptr<field>> FieldDict;