#pragma once

class Options;

template <typename BasePhysics>
class AcousticElastic2D: public BasePhysics {

 public:

  /**** Initializers ****/
  AcousticElastic2D<BasePhysics>(Options options);

};