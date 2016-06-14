#pragma once

// stl
#include <ostream>
#include <iostream>
#include <iosfwd>
#include <string>
#include <memory>

// 3rd party.
#include <mpi.h>

class Mesh;
class Options;
class Element;
class ExodusModel;

class Problem {

 public:

  virtual ~Problem() {};

  static Problem *factory(std::string solver_type);

  virtual void solve(Options options) = 0;
  virtual void initialize(Mesh *mesh, ExodusModel *model, Options &options) = 0;

};
