#pragma once

#include <Problem/Problem.h>
#include <set>

class Order2Newmark: public Problem {

  PetscReal mDt;

 public:

  /**
   * Returns individual field components for a given physical system.
   * @param physics "fluid", or "elastic2d" or "elastic3d", etc.
   * @return The defined fields (i.e. u, ux, uy, uz, ...
   */
  static std::vector<std::string> physicsToFields(const std::set<std::string> &physics);

  Order2Newmark(const std::unique_ptr<Options>& options);
  FieldDict initializeGlobalDofs(ElemVec const &elements, std::unique_ptr<Mesh> &mesh);
  FieldDict applyInverseMassMatrix(FieldDict fields);
  std::tuple<FieldDict, PetscScalar> takeTimeStep(
      FieldDict fields, PetscScalar time, std::unique_ptr<Options> const &options);

};

