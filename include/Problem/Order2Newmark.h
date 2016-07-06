#pragma once

#include <Problem/ProblemNew.h>
#include <set>

class Order2Newmark: public ProblemNew {

 public:

  FieldDict initializeGlobalDofs(ElemVec const &elements, std::unique_ptr<Mesh> &mesh);
  FieldDict applyInverseMassMatrix(FieldDict fields);
  FieldDict takeTimeStep(FieldDict fields);

};


