#include "catch.h"
#include <Eigen/Dense>
#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>
#include <stdexcept>

TEST_CASE("test functions in quadP1", "[quadP1]") {

  const static int num_dim = 2;
  const static int num_vtx = 4;

  SECTION("Interpolation on reference element") {

    Eigen::Vector4d cntr, bot_lft, bot_rgt, top_rgt, top_lft;

    cntr << 0.25, 0.25, 0.25, 0.25;
    bot_lft << 1, 0, 0, 0;
    bot_rgt << 0, 1, 0, 0;
    top_rgt << 0, 0, 1, 0;
    top_lft << 0, 0, 0, 1;

    const static int num_dim = 2;
    const static int num_vtx = 4;

    Eigen::Matrix<double, num_vtx, num_dim> vtx;
    vtx << -1, -1, 1, -1, 1, 1, -1, 1;

    REQUIRE(QuadP1::interpolateAtPoint(0, 0) == cntr);
    REQUIRE(QuadP1::interpolateAtPoint(-1, -1) == bot_lft);
    REQUIRE(QuadP1::interpolateAtPoint(1, -1) == bot_rgt);
    REQUIRE(QuadP1::interpolateAtPoint(1, 1) == top_rgt);
    REQUIRE(QuadP1::interpolateAtPoint(-1, 1) == top_lft);

  }

  SECTION("Build nodal points") {

    int num_gll = 25;
    Eigen::VectorXd pts_x_true(num_gll);
    Eigen::VectorXd pts_y_true(num_gll);

    pts_x_true
        << -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0,
        0.654654, 1, -1, -0.654654, 0, 0.654654, 1, -1, -0.654654, 0, 0.654654, 1;

    pts_y_true
        << -1, -1, -1, -1, -1, -0.654654, -0.654654, -0.654654, -0.654654, -0.654654, 0, 0, 0, 0,
        0, 0.654654, 0.654654, 0.654654, 0.654654, 0.654654, 1, 1, 1, 1, 1;


    Eigen::Matrix<double, num_vtx, num_dim> vtx;

    Eigen::VectorXd crd_r = TensorQuad<QuadP1>::GllPointsForOrder(4);
    Eigen::VectorXd crd_s = TensorQuad<QuadP1>::GllPointsForOrder(4);

    vtx << -1, -1, +1, -1, +1, +1, -1, +1;

    Eigen::VectorXd pts_x_test, pts_y_test;
    std::tie(pts_x_test, pts_y_test) = QuadP1::buildNodalPoints(crd_r, crd_s, vtx);
    REQUIRE(pts_x_test.isApprox(pts_x_true, 1e-6));
    REQUIRE(pts_y_test.isApprox(pts_y_true, 1e-6));

  }

  SECTION("Check hull") {

    Eigen::Matrix<double, num_vtx, num_dim> vtx;
    vtx << -1, -1, +1, -1, +1, +1, -1, +1;

    REQUIRE(QuadP1::checkHull(0, 0, vtx));
    REQUIRE(QuadP1::checkHull(1, 1, vtx));
    REQUIRE_FALSE(QuadP1::checkHull(-1.1, -1, vtx));

    // Make some weird element.
    vtx << -2, -2, +2, +2, +2, -2, -2, +2;
    REQUIRE_THROWS_AS(QuadP1::checkHull(0, 0, vtx), std::runtime_error);
    vtx << -1, -1, -1, +1, +1, +1, +1, -1;
    REQUIRE_THROWS_AS(QuadP1::checkHull(0, 0, vtx), std::runtime_error);

  }

  SECTION("Inverse coordinate transform") {

    Eigen::Matrix<double, num_vtx, num_dim> vtx;
    vtx << -2, -2, +2, -2, +2, +2, -2, +2;

    Eigen::Vector2d true_center, true_corner;
    true_center << 0.0, 0.0;
    true_corner << 1.0, 1.0;

    REQUIRE(QuadP1::inverseCoordinateTransform(0, 0, vtx).isApprox(true_center, 1e-6));
    REQUIRE(QuadP1::inverseCoordinateTransform(2, 2, vtx).isApprox(true_corner, 1e-6));

    // Screw up to make a weird element.
    vtx << -2, -2, +2, +2, +2, -2, -2, +2;
    REQUIRE_THROWS_AS(QuadP1::inverseCoordinateTransform(100, 100, vtx), std::runtime_error);

  }

  SECTION("Inverse Jacobian") {

    Eigen::Matrix<PetscReal, 2, 2> true_inv;
    true_inv << 0.325581, 0.0465116, 0.232558, 0.604651;
    PetscReal true_det = 5.375;

    Eigen::Matrix<PetscReal, num_vtx, num_dim> vtx;
    vtx << -1, -1, +2, -2, +3, +2, -7, +2;

    PetscReal detJac;
    Eigen::Matrix<PetscReal, 2, 2> invJac;
    std::tie(invJac, detJac) = QuadP1::inverseJacobianAtPoint(0, 0, vtx);

    REQUIRE(invJac.isApprox(true_inv, 1e-6));
    REQUIRE(detJac == Approx(true_det));

  }

}