#include "catch.h"

#include <iostream>
#include <Eigen/Dense>

#include <petsc.h>
#include <Mesh/Mesh.h>
#include <Element/Element.h>
#include <Element/ElementAdapter.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>

using namespace Eigen;
using namespace std;

template <typename Element>
class TestPlugin: public Element {

 public:
  TestPlugin<Element>(std::unique_ptr<Options> const &options): Element(options) {};

};
TEST_CASE("test_quad", "[tensor_quad]") {

  // Precision with which to compare matrices.
  double precision = 1e-5;

  // Matrix representing edge integral values.
  MatrixXd QuadEdgeTrue = MatrixXd::Zero(25, 4);

  QuadEdgeTrue(0,0)  = -235.65;
  QuadEdgeTrue(1,0)  = -204.081;
  QuadEdgeTrue(2,0)  = +28.2667;
  QuadEdgeTrue(3,0)  = -23.9857;
  QuadEdgeTrue(4,0)  = -60.15;

  QuadEdgeTrue(4,1)  = -140.35;
  QuadEdgeTrue(9,1)  = -382.667;
  QuadEdgeTrue(14,1) = -65.9556;
  QuadEdgeTrue(19,1) = +1.5562;
  QuadEdgeTrue(24,1) = +1.75;

  QuadEdgeTrue(20,2) = +1.65;
  QuadEdgeTrue(21,2) = +2.68116;
  QuadEdgeTrue(22,2) = +0.733333;
  QuadEdgeTrue(23,2) = +1.53551;
  QuadEdgeTrue(24,2) = +0.75;

  QuadEdgeTrue(0,3)  = -549.85;
  QuadEdgeTrue(5,3)  = -1481.61;
  QuadEdgeTrue(10,3) = -233.956;
  QuadEdgeTrue(15,3) = +2.89441;
  QuadEdgeTrue(20,3) = +3.85;

  // Set up custom command line arguments.
  PetscOptionsClear();
  const char *arg[] = {
      "salvus_test",
      "--testing", "true",
      "--element_shape", "quad_new",
      "--polynomial_order", "4", NULL
  };

  // Fake setting via command line.
  char **argv = const_cast<char **> (arg);
  int argc = sizeof(arg) / sizeof(const char *) - 1;
  PetscOptionsInsert(&argc, &argv, NULL);

  // Set vertices.
  Eigen::Matrix<double,4,2> vtx;
  vtx << -2, -6, +1, -6, +1, +1, -2, +1;

  SECTION("Test edge integrals") {
    // Some default set of coupling edges.
    std::vector<PetscInt> cpl = {0, 1, 2, 3};

    int max_order = 4;
    for (int order = 4; order < max_order + 1; order++) {

      std::string num = std::to_string((long long) order);
      PetscOptionsSetValue("--polynomial_order", num.c_str());

      // Initialize options.
      std::unique_ptr<Options> options(new Options);
      options->setOptions();

      // Setup test element.
      TestPlugin<TensorQuad<QuadP1>> basequad(options);
      basequad.SetVtxCrd(vtx);
      basequad.SetCplEdg(cpl);

      // Set up functions (order x**n*y**N-1)
      int ord = options->PolynomialOrder();
      Eigen::VectorXi x_exp, y_exp;
      x_exp = y_exp = Eigen::VectorXi::LinSpaced(ord + 1, 0, ord);
      y_exp(ord) = x_exp(ord - 1);

      Eigen::VectorXd pts_x, pts_y;
      std::tie(pts_x, pts_y) = QuadP1::buildNodalPoints(
          TensorQuad<QuadP1>::GllPointsForOrder(ord), TensorQuad<QuadP1>::GllPointsForOrder(ord), vtx);

      // Set up functions at GLL points.
      double x = 0, y = 0;
      Eigen::VectorXd gll_val, coords = TensorQuad<QuadP1>::GllPointsForOrder(ord);
      gll_val.setZero(basequad.NumIntPnt());
      int num_pts_p_dim = sqrt(basequad.NumIntPnt());
      for (int o = 0; o < x_exp.size(); o++) {
        for (int i = 0; i < num_pts_p_dim; i++) {
          for (int j = 0; j < num_pts_p_dim; j++) {

            int ind = i + j * num_pts_p_dim;
            gll_val(ind) += pow(pts_x(ind), x_exp(o)) * pow(pts_y(ind), y_exp(o));

          }
        }
      }

      for (int edge = 0; edge < 4; edge++) {
        REQUIRE(QuadEdgeTrue.col(edge).isApprox(basequad.applyTestAndIntegrateEdge(gll_val, edge), precision));
      }

    }
  }

}
