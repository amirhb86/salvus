#include "catch.h"
#include <Eigen/Dense>
#include <Element/HyperCube/Quad.h>
#include <Element/HyperCube/QuadP1.h>

TEST_CASE("test quadP1", "[quad/quadP1]") {

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
    vtx << -1, -1,
        1, -1,
        1, 1,
        -1, 1;

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

    pts_x_true <<
        -1, -0.654654, 0, 0.654654, 1,
        -1, -0.654654, 0, 0.654654, 1,
        -1, -0.654654, 0, 0.654654, 1,
        -1, -0.654654, 0, 0.654654, 1,
        -1, -0.654654, 0, 0.654654, 1;

    pts_y_true <<
        -1, -1, -1, -1, -1,
        -0.654654, -0.654654, -0.654654, -0.654654, -0.654654,
        0, 0, 0, 0, 0,
        0.654654, 0.654654, 0.654654, 0.654654, 0.654654,
        1, 1, 1, 1, 1;



    Eigen::Matrix<double,num_vtx,num_dim> vtx;

    Eigen::VectorXd crd_r = Quad<QuadP1>::GllPointsForOrder(4);
    Eigen::VectorXd crd_s = Quad<QuadP1>::GllPointsForOrder(4);

    vtx << -1, -1,
           +1, -1,
           +1, +1,
           -1, +1;

    Eigen::VectorXd pts_x_test, pts_y_test;
    std::tie(pts_x_test, pts_y_test) = QuadP1::buildNodalPoints(crd_r, crd_s, vtx);
    REQUIRE(pts_x_test.isApprox(pts_x_true, 1e-6));
    REQUIRE(pts_y_test.isApprox(pts_y_true, 1e-6));

  }


}