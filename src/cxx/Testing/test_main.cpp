#define CATCH_CONFIG_RUNNER
#include "catch.h"
#include <Eigen/Dense>
#include <petsc.h>
#include "../../../include/Element/ElementAdapter.h"
#include <Element/HyperCube/Quad.h>
#include <Element/HyperCube/QuadP1.h>
#include <Physics/Acoustic2D.h>
#include <Element/Simplex/Triangle.h>
#include <Element/Simplex/TriP1.h>
#include <Element/HyperCube/Quad.h>
#include <Model/ExodusModel.h>

#include <Utilities/Logging.h>

INIT_LOGGING_STATE();

int main(int argc, char *argv[]) {

  // Init Salvus command line arguments.
  PetscInitialize(&argc, &argv, NULL, NULL);

  // Run all unit tests.
  int result = Catch::Session().run(argc, argv);

  // Clean up PETSc.
  PetscFinalize();

  return result;
}

std::shared_ptr<ElementAdapter<Acoustic2D<Quad<QuadP1>>>> setup_simple_quad(Options options) {

  // Simple model.
  ExodusModel *model = new ExodusModel(options);
  model->initializeParallel();

  // Get element from options.
//  std::shared_ptr<Element> reference_element = Element::factory(options);
//  std::shared_ptr<Quad> reference_quad = std::dynamic_pointer_cast<Quad> (reference_element);
  auto reference_quad = std::make_shared<ElementAdapter<Acoustic2D<Quad<QuadP1>>>>(options);

  // Make things easy by assuming a reference element.
  // NOTE THE ELEMENT IS DISTORTED x -> [-2, 1], y -> [-6, 1]
  Eigen::Matrix<double,4,2> coord;
//  coord << -2, +1, +1, -2,
//    -6, -6, +1, +1;
  coord << -2, -6,
           +1, -6,
           +1, +1,
           -2, +1;
//  reference_quad->SetVtxCrd(coord);
//  reference_quad->attachMaterialProperties(model);
//  reference_quad->setupGradientOperator(options.PolynomialOrder());

  return reference_quad;

}

TEST_CASE("Test whether simple stuff works.", "[element]") {

  int max_order = 10;
  Eigen::VectorXd exact(max_order);
  exact <<    21/2.0, -21/4.0, -7.0, -10535/16.0,
    -1835099/400.0, -19962919/400.0,
    -177738369/400.0, -7111477851/1600.0,
    -207368760073/4800, -2094734230553/4800;

  for (int order = 1; order < max_order+1; order++) {

    std::string num = std::to_string((long long) order);
//
//    // Set up custom command line arguments.
//    PetscOptionsClear();
//    const char *arg[] = {
//        "salvus_test",
//        "--duration", "0.01",
//        "--time_step", "1e-3",
//        "--exodus_file_name", "homogeneous_iso_cartesian_2D_50s.e",
//        "--exodus_model_file_name", "homogeneous_iso_cartesian_2D_50s.e",
//        "--mesh_type", "newmark",
//        "--element_shape", "quad",
//        "--physics_system", "acoustic",
//        "--polynomial_order", num.c_str(), NULL};
//    char **argv = const_cast<char**> (arg);
//    int argc = sizeof(arg) / sizeof(const char*) - 1;
//    PetscOptionsInsert(&argc, &argv, NULL);
//    Options options;
//    options.setOptions();
//
//    std::shared_ptr<Quad<Acoustic2D<QuadP1>>> reference_quad = setup_simple_quad(options);
//
//    // Set up functions (order x**N*y**N-1)
//    int ord = options.PolynomialOrder();
//    Eigen::VectorXi x_exp = Eigen::VectorXi::LinSpaced(ord+1, 0, ord);
//    Eigen::VectorXi y_exp = x_exp;
//    y_exp[ord] = x_exp[ord - 1];
//
//    Eigen::VectorXd pts_x, pts_y;
//    std::tie(pts_x,pts_y) = QuadP1::buildNodalPoints(
//        reference_quad->IntCrdR(), reference_quad->IntCrdS(),
//        reference_quad->VtxCrd());
//
//    // Set up function values at GLL points.
//    double x, y;
//    x = y = 0.0;
//    Eigen::VectorXd gll_val;
//    Eigen::VectorXd coords = Quad::GllPointsForOrder(options.PolynomialOrder());
//    gll_val.setZero(reference_quad->NumIntPnt());
//    int num_pts_p_dim = sqrt(reference_quad->NumIntPnt());
//    for (int o = 0; o < x_exp.size(); o++) {
//      for (int i = 0; i < num_pts_p_dim; i++) {
//        for (int j = 0; j < num_pts_p_dim; j++) {
//
//          int ind = i + j * num_pts_p_dim;
//          gll_val(ind) += pow(pts_x(ind), x_exp(o)) * pow(pts_y(ind), y_exp(o));
//
//        }
//      }
//    }
//
//    // Test against analytical solution (from sympy), within floating point precision.
////    std::cout << reference_quad->integrateField(gll_val) << std::endl;
//    REQUIRE(reference_quad->integrateField(gll_val) == Approx(exact(order-1)));
  }


}

//TEST_CASE("Test triangle velocity interpolation", "[element]") {
//
//  // testing triangle with vertices (-1,-1),(1,-1),(0,sqrt(2))
//
//  Eigen::MatrixXd mMaterialVelocityAtVertices_i(3,3);
//  mMaterialVelocityAtVertices_i <<
//    1,1,2,
//    2,1,1,
//    1,2,1;
//
//  Eigen::MatrixXd check_velocity_i(3,12);
//  check_velocity_i <<
//    1.20735, 1.20735, 1.5853, 1, 1, 1.2935, 1.7065, 1.7065, 1.2935, 1, 1, 2,
//    1.5853, 1.20735, 1.20735, 1.7065,1.2935,1,1,1.2935,1.7065,2,1,1,
//    1.20735, 1.5853, 1.20735, 1.2935, 1.7065, 1.7065, 1.2935, 1, 1, 1, 2, 1;
//
//  for(int i=0;i<mMaterialVelocityAtVertices_i.rows();i++) {
//    Eigen::VectorXd mMaterialVelocityAtVertices = mMaterialVelocityAtVertices_i.row(i);
//    Eigen::VectorXd check_velocity = check_velocity_i.row(i);
//    Options options;
//    options.__SetPolynomialOrder(2);
//
//    Eigen::VectorXd mIntegrationCoordinatesR;
//    Eigen::VectorXd mIntegrationCoordinatesS;
//    std::tie(mIntegrationCoordinatesR,mIntegrationCoordinatesS) = Triangle<TriP1>::QuadraturePoints(3);
//
//    int n=0;
//    Eigen::VectorXd velocity(mIntegrationCoordinatesS.size());
//    for (auto i = 0; i < mIntegrationCoordinatesS.size(); i++) {
//
//      // Eps and eta coordinates.
//      double r = mIntegrationCoordinatesR[i];
//      double s = mIntegrationCoordinatesS[i];
//      // Get material parameters at this node.
//      auto interpolate = Triangle<TriP1>::interpolateAtPoint(r, s);
//      velocity[n] = interpolate.dot(mMaterialVelocityAtVertices);
//      n++;
//    }
//
//    REQUIRE((velocity.array()-check_velocity.array()).abs().maxCoeff() < 1e-5);
//  }
//
//}

TEST_CASE("Test quad velocity interpolation", "[element]") {

  Eigen::MatrixXd mMaterialVelocityAtVertices_i(3,4);
  mMaterialVelocityAtVertices_i <<
    1,2,2,1,
    1,1,2,2,
    1,1,2,1;
    
  Eigen::MatrixXd check_velocity_i(3,9);
  check_velocity_i <<
    1.0,1.5,2.0,1.0,1.5,2.0,1.0,1.5,2.0,
    1.0,1.0,1.0,1.5,1.5,1.5,2.0,2.0,2.0,
    1.0,1.0,1.0,1.0,1.25,1.5,1.0,1.5,2.0;
    
  for(int i=0;i<mMaterialVelocityAtVertices_i.rows();i++) {
    Eigen::VectorXd mMaterialVelocityAtVertices = mMaterialVelocityAtVertices_i.row(i);
    Eigen::VectorXd check_velocity = check_velocity_i.row(i);
    Options options;
    options.__SetPolynomialOrder(2);
    
    auto mIntegrationCoordinatesEta = Quad<QuadP1>::GllPointsForOrder(options.PolynomialOrder());
    auto mIntegrationCoordinatesEps = Quad<QuadP1>::GllPointsForOrder(options.PolynomialOrder());
    auto mNumberIntegrationPointsEta = mIntegrationCoordinatesEta.size();
    auto mNumberIntegrationPointsEps = mIntegrationCoordinatesEps.size();
    
    int n=0;
    Eigen::VectorXd velocity(mNumberIntegrationPointsEta*mNumberIntegrationPointsEps);
    for (auto eta_index = 0; eta_index < mNumberIntegrationPointsEta; eta_index++) {
      for (auto eps_index = 0; eps_index < mNumberIntegrationPointsEps; eps_index++) {

        // Eps and eta coordinates.
        double eta = mIntegrationCoordinatesEta[eta_index];
        double eps = mIntegrationCoordinatesEps[eps_index];
        // Get material parameters at this node.
        // new way
        auto interpolate1 = Quad<QuadP1>::interpolateAtPoint(eps, eta);
        velocity[n] = interpolate1.dot(mMaterialVelocityAtVertices);
        n++;
      }
    }
    
    REQUIRE((velocity.array()-check_velocity.array()).abs().maxCoeff() < 1e-5);
  }
}
