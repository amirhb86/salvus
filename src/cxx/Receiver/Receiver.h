#pragma once

#include <Utilities/Options.h>
class Receiver {

  double mPysLocX, mPysLocY, mPysLocZ;
  /** < Locations in the physical mesh */

  double mRefLocR, mRefLocS, mRefLocT;
  /** < Locations in element coordinates */

 public:

  Receiver(Options options);
  static std::vector<Receiver *> factory(Options options);

};
