#pragma once

#include <Mesh/Mesh.h>

// forward decl.

class ElasticAcousticNewmark3D: public Mesh {

 public:
  ~ElasticAcousticNewmark3D() {};
  virtual void advanceField(double dt);
  virtual void applyInverseMassMatrix();

  ElasticAcousticNewmark3D(const std::unique_ptr<Options> &options): Mesh(options) {
      mGlobalFields = {"ux",  "uy",  "uz",  "u",
                       "vx",  "vy",  "vz",  "v",
                       "ax",  "ay",  "az",  "a",
                       "ax_", "ay_", "az_", "a_",
                       "m"};}

};


