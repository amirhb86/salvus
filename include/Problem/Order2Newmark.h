#pragma once

#include <Problem/ProblemNew.h>
#include <set>

class Order2Newmark: public ProblemNew {

 public:

  Order2Newmark(const std::unique_ptr<Options>& options): ProblemNew(options) {};
  FieldDict initializeGlobalDofs(ElemVec const &elements, std::unique_ptr<Mesh> &mesh);
  FieldDict applyInverseMassMatrix(FieldDict fields);
  std::tuple<FieldDict, PetscScalar> takeTimeStep(
      FieldDict fields, PetscScalar time, std::unique_ptr<Options> const &options);

};


