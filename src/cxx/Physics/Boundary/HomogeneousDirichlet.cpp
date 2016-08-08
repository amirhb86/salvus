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
#include <Element/HyperCube/Hexahedra.h>
#include <Element/HyperCube/HexP1.h>
#include <Physics/Elastic2D.h>
#include <Physics/Elastic3D.h>
#include <Physics/AcousticElastic2D.h>
#include <Physics/ElasticAcoustic2D.h>
#include <Element/Simplex/Triangle.h>
#include <Element/Simplex/Tetrahedra.h>
#include <Element/Simplex/TriP1.h>
#include <Element/Simplex/TetP1.h>

template class HomogeneousDirichlet<
        Scalar<
            TensorQuad<
                QuadP1>>>;

template class HomogeneousDirichlet<
  Scalar<
  Triangle<
    TriP1>>>;

template class HomogeneousDirichlet<
    Elastic2D<
        TensorQuad<
            QuadP1>>>;
template class HomogeneousDirichlet<
    AcousticToElastic2D<
        Elastic2D<
            TensorQuad<
                QuadP1>>>>;
template class HomogeneousDirichlet<
    ElasticToAcoustic2D<
        Scalar<
            TensorQuad<
                QuadP1>>>>;
template class HomogeneousDirichlet<
    Scalar<
        Hexahedra<
            HexP1>>>;
template class HomogeneousDirichlet<
    Elastic3D<
        Hexahedra<
            HexP1>>>;
