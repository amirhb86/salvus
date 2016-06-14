#pragma once

// stl.
#include <vector>
#include <memory>

// forward decl.
class Options;

class Source {

  double mPhysicalLocationX;
  double mPhysicalLocationY;
  double mPhysicalLocationZ;

  double mReferenceLocationR;
  double mReferenceLocationS;
  double mReferenceLocationT;

 public:

  static std::vector<std::shared_ptr<Source>> factory(std::unique_ptr<Options> const &options);

  /**
   * Set the source x coordinate in physical (model) space.
   */
  inline void SetPhysicalLocationX(double location_x) { mPhysicalLocationX = location_x; }

  /**
   * Set the source y coordinate in physical (model) space.
   */
  inline void SetPhysicalLocationY(double location_y) { mPhysicalLocationY = location_y; }

  /**
   * Set the source z coordinate in physical (model) space.
   */
  inline void SetPhysicalLocationZ(double location_z) { mPhysicalLocationZ = location_z; }

  /**
   * Set the source epsilon coordinate on the reference element.
   */
  inline void setReferenceLocationR(double location_r) { mReferenceLocationR = location_r; }

  /**
   * Set the source MAYBE WE SHOULD CHANGE COORDINATE NAMES coordinate on the reference element.
   */
  inline void setReferenceLocationS(double location_s) { mReferenceLocationS = location_s; }

  /**
   * Set the source eta coordinate on the reference element.
   */
  inline void setReferenceLocationT(double location_t) { mReferenceLocationT = location_t; }

  inline double PhysicalLocationX() { return mPhysicalLocationX; }
  inline double PhysicalLocationY() { return mPhysicalLocationY; }
  inline double PhysicalLocationZ() { return mPhysicalLocationZ; }

  inline double ReferenceLocationR() { return mReferenceLocationR; }
  inline double ReferenceLocationS() { return mReferenceLocationS; }
  inline double ReferenceLocationT() { return mReferenceLocationT; }

  /**
   * Returns a value for the force, given a certain time. This needs to be implemented by each derived class. For
   * instance, a Ricker source will need to implement the source time characeristics of a ricker source time
   * function.
   */
  virtual double fire(const double &time) = 0;

};

class Ricker: public Source {

  double mAmplitude;
  double mTimeDelay;
  double mCenterFreq;

 public:

  Ricker(std::unique_ptr<Options> const &options, int number);
  double fire(const double &time);

};

