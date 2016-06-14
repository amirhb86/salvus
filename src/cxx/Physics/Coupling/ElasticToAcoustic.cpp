#include <Model/ExodusModel.h>
#include <Mesh/Mesh.h>
#include <Utilities/Options.h>
#include <Physics/ElasticAcoustic2D.h>

using namespace Eigen;

template <typename BasePhysics>
ElasticToAcoustic2D<BasePhysics>::ElasticToAcoustic2D(std::unique_ptr<Options> const &options): BasePhysics(options) { }

template <typename BasePhysics>
std::vector<std::string> ElasticToAcoustic2D<BasePhysics>::PullElementalFields() const {
  return {"u", "vx", "vy"};
}

template <typename BasePhysics>
void ElasticToAcoustic2D<BasePhysics>::setBoundaryConditions(Mesh *mesh) {
  for (auto e: mesh->CouplingFields(BasePhysics::ElmNum())) {
    mEdg.push_back(std::get<0>(e));
    mNbr.push_back(mesh->GetNeighbouringElement(mEdg.back(), BasePhysics::ElmNum()));
    mNbrCtr.push_back(mesh->getElementCoordinateClosure(mNbr.back()).colwise().mean());
  }
  BasePhysics::setBoundaryConditions(mesh);
}

template <typename BasePhysics>
MatrixXd ElasticToAcoustic2D<BasePhysics>::computeSurfaceIntegral(const Ref<const MatrixXd> &u) {

  // col0->potential, col1->ux, col2->uy.
  MatrixXd rval = Eigen::MatrixXd::Zero(BasePhysics::NumIntPnt(), 1);

  Vector2d n;
  for (auto e: mEdg) {
    n = BasePhysics::getEdgeNormal(e);
    VectorXd transform = u.rightCols(2) * n;
    rval += BasePhysics::applyTestAndIntegrateEdge(transform, e);
  }

  return rval;

}

#include <Physics/Acoustic2D.h>
#include <Element/HyperCube/Quad.h>
#include <Element/HyperCube/QuadP1.h>
template class ElasticToAcoustic2D<Acoustic2D<Quad<QuadP1>>>;
