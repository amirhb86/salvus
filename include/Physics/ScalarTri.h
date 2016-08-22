#pragma once

// stl.
#include <vector>

// 3rd party.
#include <Eigen/Dense>
#include <Utilities/Types.h>

// forward decl.
class Mesh;
class Options;
class ExodusModel;

template <typename Shape>
class ScalarTri: public Shape {
  /**
   * \class ScalarTri
   *
   * Custom class for triangles
   */

 private:

  RealVec mStiff;

 public:

  /**** Initializers ****/
  ScalarTri<Shape>(std::unique_ptr<Options> const &options): Shape(options) {};
  ~ScalarTri<Shape>() {};
  
  /**** Time loop functions ****/
  
  RealMat computeStiffnessTerm(const Eigen::Ref<const RealMat>& u);

  const static std::string Name() { return "ScalarTri_" + Shape::Name(); }

};

