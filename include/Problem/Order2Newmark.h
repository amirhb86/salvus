#pragma once

#include <Problem/Problem.h>
#include <set>

class Order2Newmark: public Problem {

  PetscReal mDt;

 public:

  Order2Newmark(const std::unique_ptr<Options>& options);
  FieldDict initializeGlobalDofs(ElemVec const &elements, std::unique_ptr<Mesh> &mesh);
  FieldDict applyInverseMassMatrix(FieldDict fields);
  std::tuple<FieldDict, PetscScalar> takeTimeStep(
      FieldDict fields, PetscScalar time, std::unique_ptr<Options> const &options);

};

