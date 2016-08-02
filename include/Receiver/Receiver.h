#pragma once

// stl.
#include <memory>
#include <map>

// 3rd party.
#include <Eigen/Dense>

// forward decl.
class Options;

class Receiver {

  static PetscInt mNumRecs;
  /** < Counter for number of receivers. */

  static std::vector<std::string> names;
  /** < Static vector holding all receiver names. */

  long mNum;
  /** < Number of this receiver. */

  double mLocX, mLocY, mLocZ;
  /** < Locations in the physical mesh */

  double mRefLocR, mRefLocS, mRefLocT;
  /** < Locations in element coordinates */

  std::string mName;
  /** < Receiver name */

 protected:

  std::map<std::string,std::vector<float>> store;
  /** < Store all fields in here */

 public:

  inline void SetNum(const PetscInt num) { mNum = num; }
  Receiver(std::unique_ptr<Options> const &options);
  virtual ~Receiver();
  static std::vector<std::unique_ptr<Receiver>> Factory(std::unique_ptr<Options> const &options);

  /* Get number of active receivers. */
  static PetscInt NumReceivers() { return mNumRecs; }

  inline void SetRefLocR (double val) { mRefLocR = val; }
  inline void SetRefLocS (double val) { mRefLocS = val; }
  inline void SetRefLocT (double val) { mRefLocT = val; }

  inline long Num() { return mNum; }
  inline std::string Name() { return mName; }
  inline PetscReal LocX() { return mLocX; }
  inline PetscReal LocY() { return mLocY; }
  inline PetscReal LocZ() { return mLocZ; }
  inline PetscReal RefLocR() { return mRefLocR; }
  inline PetscReal RefLocS() { return mRefLocS; }
  inline PetscReal RefLocT() { return mRefLocT; }

  void record(const double val, const std::string &field);

  virtual void write() = 0;

};
