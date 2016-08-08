#pragma once

// stl.
#include <iostream>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>

#include <Utilities/Types.h>

class HexP1 {

 private:

  static inline PetscReal dn0dr(const PetscReal r) { return (-1) * (1 - r) / 4.0; }
  static inline PetscReal dn1dr(const PetscReal r) { return (+1) * (1 - r) / 4.0; }
  static inline PetscReal dn2dr(const PetscReal r) { return (+1) * (1 + r) / 4.0; }
  static inline PetscReal dn3dr(const PetscReal r) { return (-1) * (1 + r) / 4.0; }
  static inline PetscReal dn0ds(const PetscReal s) { return (1 - s) * -1.0 / 4.0; }
  static inline PetscReal dn1ds(const PetscReal s) { return (1 + s) * -1.0 / 4.0; }
  static inline PetscReal dn2ds(const PetscReal s) { return (1 + s) * 1.0 / 4.0; }
  static inline PetscReal dn3ds(const PetscReal s) { return (1 - s) * 1.0 / 4.0; }

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
  static void faceJacobianAtPoint(
      const PetscReal r, const PetscReal s, const Eigen::Ref<const QuadVtx> &vtx,
      PetscReal &detJac);

  static std::tuple<RealVec, RealVec, RealVec>
  buildNodalPoints(const Eigen::Ref<const RealVec>& intCrdR,
                   const Eigen::Ref<const RealVec>& intCrdS,
                   const Eigen::Ref<const RealVec>& intCrdT,
                   const Eigen::Ref<const HexVtx>& vtx);

  const static std::string Name() { return "HexP1"; }

  /// Class name
  const static ElementType type() { return ElementType::HEXP1; }
};


