#include <Mesh/Mesh.h>
#include <Physics/HomogeneousDirichlet.h>

using namespace Eigen;

template <typename Base>
HomogeneousDirichlet<Base>::HomogeneousDirichlet(
    std::unique_ptr<Options> const &options): Base(options) {

}

template <typename Base>
void HomogeneousDirichlet<Base>::setBoundaryConditions(std::unique_ptr<Mesh> const &mesh) {

  auto bnds = mesh->BoundaryPoints();
  PetscInt num_pts; DMPlexGetConeSize(mesh->DistributedMesh(), Base::ElmNum(), &num_pts);
  const PetscInt *pts = NULL; DMPlexGetCone(mesh->DistributedMesh(), Base::ElmNum(), &pts);
  for (PetscInt i = 0; i < num_pts; i++) {
    if (bnds.find(pts[i]) != bnds.end()) { mBndEdg.push_back(i); }
  }

  Base::setBoundaryConditions(mesh);

}

template <typename Base>
RealMat HomogeneousDirichlet<Base>::computeStiffnessTerm(const Ref<const RealMat> &u) {

  RealMat s = Base::computeStiffnessTerm(u);
  for (PetscInt i = 0; i < s.cols(); i++) {
    for (auto edge: mBndEdg) { Base::setEdgeToValue(edge, 0, s.col(i)); }
  }

  return s;

}



#include <Physics/Scalar.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>
template class HomogeneousDirichlet<
        Scalar<
            TensorQuad<
                QuadP1>>>;
