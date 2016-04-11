#pragma once

#include <Utilities/Options.h>
#include <memory>
class Receiver {

  static long num;
  /** < Counter for number of receivers. */

  long mNum;
  /** < Number of this receiver. */

  double mPysLocX1, mPysLocX2, mPysLocX3;
  /** < Locations in the physical mesh */

  double mRefLocR, mRefLocS, mRefLocT;
  /** < Locations in element coordinates */

 public:

  // See discussion I've had with myself in the Element header (mRec).
  virtual std::shared_ptr<Receiver> clone() const = 0;

  Receiver(Options options);
  virtual ~Receiver();
  static std::vector<std::shared_ptr<Receiver>> factory(Options options);

  inline void SetRefLocR (double val) { mRefLocR = val; }
  inline void SetRefLocS (double val) { mRefLocS = val; }
  inline void SetRefLocT (double val) { mRefLocT = val; }

  inline long Num() { return mNum; }
  inline double PysLocX1() { return mPysLocX1; }
  inline double PysLocX2() { return mPysLocX2; }
  inline double PysLocX3() { return mPysLocX3; }

};
