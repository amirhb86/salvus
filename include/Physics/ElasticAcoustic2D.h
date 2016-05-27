#pragma once

class Options;

template <typename BasePhysics>
class ElasticAcoustic2D: public BasePhysics {

 public:

  /**** Initializers ****/
  ElasticAcoustic2D<BasePhysics>(Options options);

};
