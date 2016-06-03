#include <Physics/AcousticElastic2D.h>
#include <Utilities/Options.h>
#include <Mesh/Mesh.h>

using namespace Eigen;

template <typename BasePhysics>
AcousticToElastic2D<BasePhysics>::AcousticToElastic2D(Options options): BasePhysics(options) {

//  uElastic.setZero(BasePhysics::NumIntPnt(), num_elastic_components);
//  uScalar.setZero(BasePhysics::NumIntPnt(), num_scalar_components);

}

template <typename BasePhysics>
std::vector<std::string> AcousticToElastic2D<BasePhysics>::PullElementalFields() const {
  return {"ux", "uy", "v"};
}

template <typename BasePhysics>
void AcousticToElastic2D<BasePhysics>::setBoundaryConditions(Mesh *mesh) {
  for (auto e: mesh->CouplingFields(BasePhysics::ElmNum())) {
    mEdg.push_back(std::get<0>(e));
  }
  BasePhysics::setBoundaryConditions(mesh);
}

template <typename BasePhysics>
Eigen::MatrixXd AcousticToElastic2D<BasePhysics>::computeSurfaceIntegral(const Eigen::Ref<const Eigen::MatrixXd> &u) {

  // col0->ux, col1->uy, col2->potential.
  Eigen::VectorXd rval = Eigen::VectorXd::Zero(BasePhysics::NumIntPnt());
  for (auto e : mEdg){
    rval.noalias() += BasePhysics::applyTestAndIntegrateEdge(u.col(2), e);
  }
  return 1000 * rval;

}

#include <Physics/Acoustic2D.h>
#include <Element/HyperCube/Quad.h>
#include <Element/HyperCube/QuadP1.h>
template class AcousticToElastic2D<Acoustic2D<Quad<QuadP1>>>;