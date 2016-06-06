#pragma once

#include <Mesh/Mesh.h>

class ScalarLeapFrog4th : public Mesh {

public:


  ~ScalarLeapFrog4th() {};
    
  ScalarLeapFrog4th() { mGlobalFields = {"u", "unm1", "unp1", "a0", "a1", "m", "u_exact"}; }
  ScalarLeapFrog4th(std::vector<std::string> fields) {mGlobalFields = fields;}
  
  virtual void advanceField(double dt);
  virtual void applyInverseMassMatrix();

  void applyInverseMassMatrix(std::string tofield);
  
};
