#include "catch.h"
#include <Eigen/Dense>
#include <Utilities/Types.h>
#include <Mesh/Mesh.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Physics/Scalar.h>
#include <Element/ElementAdapter.h>
#include <Problem/Problem.h>
#include <petscviewerhdf5.h>
#include <Element/HyperCube/Hexahedra.h>
#include <Element/HyperCube/HexP1.h>

#include <stdexcept>

TEST_CASE("test functions in HexP1", "[hexP1]") {

  const static PetscInt num_dim = 3;
  const static PetscInt num_vtx = 8;

  SECTION("Interpolation on reference element") {

    Eigen::Matrix<PetscReal,8,1> ctr;
    Eigen::Matrix<PetscReal,8,1> bot_lft_bk, bot_rgt_bk, bot_lft_ft, bot_rgt_ft;
    Eigen::Matrix<PetscReal,8,1> top_lft_bk, top_rgt_bk, top_lft_ft, top_rgt_ft;

    ctr << 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125;
    bot_lft_bk << 1, 0, 0, 0, 0, 0, 0, 0;
    bot_lft_ft << 0, 1, 0, 0, 0, 0, 0, 0;
    bot_rgt_ft << 0, 0, 1, 0, 0, 0, 0, 0;
    bot_rgt_bk << 0, 0, 0, 1, 0, 0, 0, 0;
    top_lft_bk << 0, 0, 0, 0, 1, 0, 0, 0;
    top_rgt_bk << 0, 0, 0, 0, 0, 1, 0, 0;
    top_rgt_ft << 0, 0, 0, 0, 0, 0, 1, 0;
    top_lft_ft << 0, 0, 0, 0, 0, 0, 0, 1;

    /* Inward normals. */
    REQUIRE(HexP1::interpolateAtPoint(0, 0, 0) == ctr);
    REQUIRE(HexP1::interpolateAtPoint(-1, -1, -1) == bot_lft_bk);
    REQUIRE(HexP1::interpolateAtPoint(-1, +1, -1) == bot_lft_ft);
    REQUIRE(HexP1::interpolateAtPoint(+1, +1, -1) == bot_rgt_ft);
    REQUIRE(HexP1::interpolateAtPoint(+1, -1, -1) == bot_rgt_bk);
    REQUIRE(HexP1::interpolateAtPoint(-1, -1, +1) == top_lft_bk);
    REQUIRE(HexP1::interpolateAtPoint(+1, -1, +1) == top_rgt_bk);
    REQUIRE(HexP1::interpolateAtPoint(+1, +1, +1) == top_rgt_ft);
    REQUIRE(HexP1::interpolateAtPoint(-1, +1, +1) == top_lft_ft);

  }

  SECTION("Build nodal points") {

    PetscInt num_gll = 125;
    RealVec pts_x_true(num_gll);
    RealVec pts_y_true(num_gll);
    RealVec pts_z_true(num_gll);

    pts_x_true
        << -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0,
        0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1,
        -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1, -0.654654,
        0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1,
        -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0,
        0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1,
        -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0,
        0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1,
        -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654,
        1, -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1;

    pts_y_true
        << -1, -1, -1, -1, -1, -0.654654, -0.654654, -0.654654, -0.654654, -0.654654, 0, 0, 0,
        0, 0, 0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 1, 1, 1, 1, 1, -1, -1, -1,
        -1, -1, -0.654654, -0.654654, -0.654654, -0.654654, -0.654654, 0, 0, 0, 0, 0,
        0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 1, 1, 1, 1, 1, -1, -1, -1, -1,
        -1, -0.654654, -0.654654, -0.654654, -0.654654, -0.654654, 0, 0, 0, 0, 0, 0.654654,
        0.654654, 0.654654, 0.654654, 0.654654, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1,
        -0.654654, -0.654654, -0.654654, -0.654654, -0.654654, 0, 0, 0, 0, 0, 0.654654,
        0.654654, 0.654654, 0.654654, 0.654654, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -0.654654,
        -0.654654, -0.654654, -0.654654, -0.654654, 0, 0, 0, 0, 0, 0.654654, 0.654654, 0.654654,
        0.654654, 0.654654, 1, 1, 1, 1, 1;

    pts_z_true
        << -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -0.654654, -0.654654, -0.654654, -0.654654, -0.654654,
        -0.654654, -0.654654, -0.654654, -0.654654, -0.654654, -0.654654, -0.654654,
        -0.654654, -0.654654, -0.654654, -0.654654, -0.654654, -0.654654, -0.654654,
        -0.654654, -0.654654, -0.654654, -0.654654, -0.654654, -0.654654, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.654654,
        0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 0.654654,
        0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 0.654654,
        0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;

    RealVec crd_r(5); crd_r << -1, -0.654654, 0, 0.654654, 1;
    RealVec crd_s(5); crd_s << -1, -0.654654, 0, 0.654654, 1;
    RealVec crd_t(5); crd_t << -1, -0.654654, 0, 0.654654, 1;

    HexVtx vtx;
    vtx <<
        -1, -1, -1,
        -1, +1, -1,
        +1, +1, -1,
        +1, -1, -1,
        -1, -1, +1,
        +1, -1, +1,
        +1, +1, +1,
        -1, +1, +1;

    RealVec pts_x_test, pts_y_test, pts_z_test;
    std::tie(pts_x_test, pts_y_test, pts_z_test) = HexP1::buildNodalPoints(
        crd_r, crd_s, crd_t, vtx);
    REQUIRE(pts_x_test.isApprox(pts_x_true, 1e-6));
    REQUIRE(pts_y_test.isApprox(pts_y_true, 1e-6));
    REQUIRE(pts_z_test.isApprox(pts_z_true, 1e-6));

  }

  SECTION("Check hull") {

    HexVtx vtx;
    vtx <<
        -1, -1, -1,
        -1, +1, -1,
        +1, +1, -1,
        +1, -1, -1,
        -1, -1, +1,
        +1, -1, +1,
        +1, +1, +1,
        -1, +1, +1;

    REQUIRE(HexP1::checkHull(0, 0, 0, vtx));
    REQUIRE(HexP1::checkHull(1, 1, 1, vtx));
    REQUIRE_FALSE(HexP1::checkHull(-1.1, -1, -1, vtx));

    // Make some weird element.
    vtx <<
        -1, -1, -1,
        +1, -1, +1,
        -1, +1, -1,
        +1, +1, -1,
        +1, -1, -1,
        -1, -1, +1,
        +1, +1, +1,
        -1, +1, +1;
    REQUIRE_THROWS_AS(HexP1::checkHull(0, 0, 0, vtx), std::runtime_error);

  }

  SECTION("Inverse coordinate transform") {

    HexVtx vtx;
    vtx <<
        -2, -2, -2,
        -2, +2, -2,
        +2, +2, -2,
        +2, -2, -2,
        -2, -2, +2,
        +2, -2, +2,
        +2, +2, +2,
        -2, +2, +2;

    RealVec3 true_center, true_corner;
    true_center << 0.0, 0.0, 0.0;
    true_corner << 1.0, 1.0, 1.0;

    REQUIRE(HexP1::inverseCoordinateTransform(0, 0, 0, vtx).isApprox(true_center, 1e-6));
    REQUIRE(HexP1::inverseCoordinateTransform(2, 2, 2, vtx).isApprox(true_corner, 1e-6));

    // Screw up to make a weird element.
    vtx <<
        -2, -3, -7,
        +2, -3, +7,
        -2, +3, -7,
        +2, +3, -7,
        +2, -3, -7,
        -2, -3, +7,
        +2, +3, +7,
        -2, +3, +7;
    REQUIRE_THROWS_AS(HexP1::inverseCoordinateTransform(100, 100, 100, vtx), std::runtime_error);

  }

  SECTION("Inverse Jacobian") {

    RealMat3x3 true_inv;
    true_inv <<
        0.5, -0.05, 0.05,
        0.0, +0.45, 0.05,
        0.0, +0.05, 0.45;
    PetscReal true_det = 10;

    Eigen::Matrix<PetscReal, num_vtx, num_dim> vtx;
    vtx <<
        -2, -2, -2,
        -3, +2, -2,
        +3, +3, -4,
        +2, -2, -2,
        -1, -3, +2,
        +1, -2, +2,
        +2, +2, +2,
        -2, +2, +2;

    PetscReal detJac;
    RealMat3x3 invJac;
    HexP1::inverseJacobianAtPoint(0, 0, 0, vtx, detJac, invJac);

    REQUIRE(invJac.isApprox(true_inv, 1e-6));
    REQUIRE(detJac == Approx(true_det));

  }

}
