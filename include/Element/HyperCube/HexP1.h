#pragma once

// stl.
#include <iostream>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>

#include <Utilities/Types.h>

class HexP1 {

 private:

  const static PetscInt mNumDim = 3;
  const static PetscInt mNumVtx = 8;

 public:

  HexP1() {};
  static RealVec interpolateAtPoint(
      const double r, const PetscReal s, const PetscReal t);
  static bool checkHull(
      const PetscReal x, const PetscReal y, const PetscReal z,
      const Eigen::Ref<const HexVtx> &vtx);
  static RealVec3 inverseCoordinateTransform(
      const PetscReal x, const PetscReal y, const PetscReal z,
      const Eigen::Ref<const HexVtx> &vtx);
  static void inverseJacobianAtPoint(
      const PetscReal r, const PetscReal s, const PetscReal t,
      const Eigen::Ref<const HexVtx> &vtx, PetscReal &detJac,
      Eigen::Ref<RealMat3x3> invJac);
  
  static std::tuple<RealVec, RealVec, RealVec>
  buildNodalPoints(const Eigen::Ref<const RealVec>& intCrdR,
                   const Eigen::Ref<const RealVec>& intCrdS,
                   const Eigen::Ref<const RealVec>& intCrdT,
                   const Eigen::Ref<const HexVtx>& vtx);
  
};


