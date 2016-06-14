#pragma once

// stl.
#include <vector>
#include <iostream>
#include <memory>

// salvus.
#include <Element/Element.h>

class Mesh;
class Model;
class Options;

class ProblemNew {

 public:

  ProblemNew() {};
  virtual ~ProblemNew() {};
  static std::unique_ptr<ProblemNew> Factory(std::unique_ptr<Options> const &options);

  std::vector<std::unique_ptr<Element>> initializeElements(std::unique_ptr<Mesh> &mesh,
                                                           std::unique_ptr<Model> &model,
                                                           std::unique_ptr<Options> const &options);
  void integrateElements(std::unique_ptr<Mesh> mesh);

  /* Depends on time stepping scheme. */
  virtual void applyInverseMassMatrix() = 0;
  virtual void takeTimeStep() = 0;

};


