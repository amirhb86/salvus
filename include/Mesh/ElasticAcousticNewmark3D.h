#pragma once

#include <Mesh/Mesh.h>

// forward decl.

class ElasticAcousticNewmark2D: public Mesh {

 public:
  ~ElasticAcousticNewmark2D() {};
  virtual void advanceField(double dt);
  virtual void applyInverseMassMatrix();

  ElasticAcousticNewmark2D() { mGlobalFields = {"ux",  "uy",  "uz",  "u",
                                                "vx",  "vy",  "vz",  "v",
                                                "ax",  "ay",  "az",  "a",
                                                "ax_", "ay_", "az_", "a_",
                                                "m"};}

};


