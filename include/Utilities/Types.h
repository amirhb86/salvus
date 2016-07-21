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
/// Floating points.
typedef Eigen::Matrix<PetscReal, 2, 1> RealVec2;
typedef Eigen::Matrix<PetscReal, 3, 1> RealVec3;
typedef Eigen::Matrix<PetscReal, 4, 1> RealVec4;
typedef Eigen::Matrix<PetscReal, 4, 2> QuadVtx;
typedef Eigen::Matrix<PetscReal, 8, 3> HexVtx;
typedef Eigen::Matrix<PetscReal, 2, 2> RealMat2x2;
typedef Eigen::Matrix<PetscReal, 3, 3> RealMat3x3;
typedef Eigen::Matrix<PetscReal, Eigen::Dynamic, 1> RealVec;
typedef Eigen::Matrix<PetscReal, Eigen::Dynamic, Eigen::Dynamic> RealMat;

/// Integers and idicies.
typedef Eigen::Matrix<PetscInt, Eigen::Dynamic, 1> IntVec;

/// Complex objects.
typedef std::vector<std::unique_ptr<Element>> ElemVec;
typedef std::map<std::string, std::unique_ptr<field>> FieldDict;

/// Strongly typed Element types
enum class ElementType {QUADP1, TRIP1, HEXP1, TETP1};
