#pragma once

// stl.
#include <vector>
#include <memory>

// 3rd party.
#include <petsc.h>

// forward decl.
class Options;

class Source {

  double mLocX;
  double mLocY;
  double mLocZ;

  double mLocR;
  double mLocS;
  double mLocT;

  PetscInt mNum;
  static PetscInt number;

 public:

  /* Get number of active sources. */
  static PetscInt NumSources() { return number; }

  /* Constructor. */
  Source(std::unique_ptr<Options> const &options);
  virtual ~Source();

  /* Factory method returning single sources. */
  static std::vector<std::unique_ptr<Source>> Factory(std::unique_ptr<Options> const &options);

  /* Unique numbers. */
  inline void SetNum(const PetscInt num) { mNum = num; }
  PetscInt Num() { return mNum; }

  /* Physical coordinates. */
  inline void SetLocX(double location_x) { mLocX = location_x; }
  inline void SetLocY(double location_y) { mLocY = location_y; }
  inline void SetLocZ(double location_z) { mLocZ = location_z; }

  /* Reference coordinates. */
  inline void SetLocR(double location_r) { mLocR = location_r; }
  inline void SetLocS(double location_s) { mLocS = location_s; }
  inline void SetLocT(double location_t) { mLocT = location_t; }

  inline double LocX() { return mLocX; }
  inline double LocY() { return mLocY; }
  inline double LocZ() { return mLocZ; }

  inline double LocR() { return mLocR; }
  inline double LocS() { return mLocS; }
  inline double LocT() { return mLocT; }

  /**
   * Returns a value for the force, given a certain time. This needs to be implemented by each derived class. For
   * instance, a Ricker source will need to implement the source time characteristics of a ricker source time
   * function.
   */
  virtual double fire(const double &time) = 0;

};
