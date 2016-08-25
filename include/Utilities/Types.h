#pragma once

#include <hdf5.h>
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

/// Custom exceptions.
class salvus_warning: public std::exception {
  std::string mMsg;
 public:
  salvus_warning(const std::string &msg) {
    mMsg = msg;
  }
  virtual const char* what() const throw() { return mMsg.c_str(); }
};

/// Wavefield container
template <typename T>
class UncompressedWavefieldContainer {
  PetscInt nElm, nPnt, nTsp, nCmp;
  std::vector<T> _data;
 public:
  UncompressedWavefieldContainer(PetscInt nTsp=0, PetscInt nElm=0,
                                 PetscInt nCmp=0, PetscInt nPnt=0) :
      nTsp(nTsp), nElm(nElm), nCmp(nCmp), nPnt(nPnt), _data(nTsp*nElm*nCmp*nPnt) {}
  T &operator()(PetscInt tsp, PetscInt elm, PetscInt cmp, PetscInt pnt) {
    return _data[tsp*nElm*nCmp*nPnt + elm*nCmp*nPnt + cmp*nPnt + pnt];
  }
  T &operator()(PetscInt i) {
    return _data[i];
  }
  T &data() { return _data[0]; }
  void resize(PetscInt nTsp_, PetscInt nElm_, PetscInt nCmp_, PetscInt nPnt_) {
    nElm = nElm_; nTsp = nTsp_; nCmp = nCmp_; nPnt = nPnt_;
    _data.resize(nElm*nTsp*nCmp*nPnt);
  };
  PetscInt size() { return nElm*nTsp*nCmp*nPnt; }
  PetscInt elm() { return nElm; }
  PetscInt tsp() { return nTsp; }
  PetscInt cmp() { return nCmp; }
  PetscInt pnt() { return nPnt; }
  static hid_t Hdf5Datatype();
};

template<>
inline hid_t UncompressedWavefieldContainer<double>::Hdf5Datatype() { return H5T_NATIVE_DOUBLE; }




/// Strongly typed Element types
enum class ElementType {QUADP1, TRIP1, HEXP1, TETP1};
