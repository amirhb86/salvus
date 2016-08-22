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

  PetscInt num, *pts = NULL;
  DMPlexGetTransitiveClosure(mesh->DistributedMesh(), Base::ElmNum(), PETSC_TRUE,
                             &num, &pts);
  /* Go through all graph depths (vertex -> element[d-1]) */
  for (PetscInt dep = 0, dep_c = 0; dep < Base::NumDim(); dep++, dep_c = 0) {
    /* Go through the graph closure for this element. */
    for (PetscInt i = 0; i < 2 * num; i += 2) {
      /* Get indices specifying points referring to current depth. */
      PetscInt p_beg, p_end;
      DMPlexGetDepthStratum(mesh->DistributedMesh(), dep, &p_beg, &p_end);
      /* If we're outside this depth, just continue. */
      if (! (pts[i] >= p_beg && pts[i] < p_end)) continue;
      /* If we're inside, look through all boundary points. */
      for (auto &pnt: bnds) {
        /* If the point belonging to this element is on a mesh boundary. */
        if (std::get<1>(pnt) == pts[i]) {
          /* Save the proper dof. */
          switch (dep) {
            case 0: /* vertex */
            {
              auto pv = Base::getDofsOnVtx(dep_c);
              mBndDofs.push_back(pv);
              break;
            }
            case 1: /* edge */
            {
              auto pe = Base::getDofsOnEdge(dep_c);
              mBndDofs.insert(mBndDofs.end(), pe.begin(), pe.end());
              break;
            }
            case 2: /* face */
            {
              auto pf = Base::getDofsOnFace(dep_c);
              mBndDofs.insert(mBndDofs.end(), pf.begin(), pf.end());
              break;
            }
            default:
              break;
          }
        }
      }
      dep_c++;
    }
  }
  DMPlexRestoreTransitiveClosure(mesh->DistributedMesh(), Base::ElmNum(), PETSC_TRUE,
                                 &num, &pts);

  /* Only need unique boundary points. */
  mBndDofs.erase(std::unique(mBndDofs.begin(), mBndDofs.end()), mBndDofs.end());
  Base::setBoundaryConditions(mesh);
}

template <typename Base>
RealMat HomogeneousDirichlet<Base>::computeStiffnessTerm(const Ref<const RealMat> &u) {

  RealMat s = Base::computeStiffnessTerm(u);
  for (PetscInt i = 0; i < s.cols(); i++) {
    for (auto dof: mBndDofs) { s(dof,i) = 0.0; }
  }
  return s;
}



#include <Physics/Scalar.h>
#include <Physics/ScalarTri.h>
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
  ScalarTri<
  Scalar<
  Triangle<
  TriP1>>>>;

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
