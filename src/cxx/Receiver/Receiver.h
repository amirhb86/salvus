#pragma once

#include <Utilities/Options.h>
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

  Receiver(Options options);
  ~Receiver();
  static std::vector<std::unique_ptr<Receiver>> factory(Options options);

  inline long Num() { return mNum; }
  inline double PysLocX1() { return mPysLocX1; }
  inline double PysLocX2() { return mPysLocX2; }
  inline double PysLocX3() { return mPysLocX3; }

};
