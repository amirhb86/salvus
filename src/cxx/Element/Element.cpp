#include <stdexcept>

#include <Element/Element.h>
#include <Element/ElementAdapter.h>

#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>
#include <Element/Simplex/Triangle.h>
#include <Element/Simplex/TriP1.h>
#include <Element/Simplex/Tetrahedra.h>
#include <Element/Simplex/TetP1.h>
#include <Element/HyperCube/Hexahedra.h>
#include <Element/HyperCube/HexP1.h>

#include <Physics/Scalar.h>
#include <Physics/Elastic2D.h>
#include <Physics/Elastic3D.h>
#include <Physics/AcousticElastic2D.h>
#include <Physics/ElasticAcoustic2D.h>
#include <Utilities/Utilities.h>
#include <Utilities/Logging.h>
#include <Physics/HomogeneousDirichlet.h>

/* Define all possible element classes as types here. */
typedef class ElementAdapter<Scalar<TensorQuad<QuadP1>>> ScalarQuadP1;
typedef class ElementAdapter<Elastic2D<TensorQuad<QuadP1>>> ElasticQuadP1;

/* Coupled classes. */
typedef class ElementAdapter<AcousticToElastic2D<Elastic2D<TensorQuad<QuadP1>>>> AcousticCplElasticQuadP1;
typedef class ElementAdapter<ElasticToAcoustic2D<Scalar<TensorQuad<QuadP1>>>> ElasticCplAcousticQuadP1;

/* Boundary conditions. */
typedef class ElementAdapter<HomogeneousDirichlet<Scalar<TensorQuad<QuadP1>>>> Tester;

std::unique_ptr<Element> Element::Factory(
    const std::vector<std::string>& physics_base,
    const std::vector<std::string>& physics_couple,
    std::unique_ptr<Options> const &options) {

  // define field combinations.
  std::vector<std::string> fluid {"fluid"};
  std::vector<std::string> elastic_2d {"2delastic"};
  std::vector<std::string> elastic_3d {"3delastic"};
  std::vector<std::string> boundary {"boundary"};

  // find correct combination of coupling physics.
  std::set<std::string> coupling_set(physics_couple.begin(), physics_couple.end());
  bool fl, el2d, el3d, bnd;
  fl    = coupling_set.find("fluid")      != coupling_set.end() ? true : false;
  el2d  = coupling_set.find("2delastic")  != coupling_set.end() ? true : false;
  el3d  = coupling_set.find("3delastic")  != coupling_set.end() ? true : false;
  bnd   = coupling_set.find("boundary")   != coupling_set.end() ? true : false;

  try {

    if (options->ElementShape() == "quad_new") {
      /* If only acoustic, return a base acoustic. */
      if (physics_base == fluid) {
        if (!physics_couple.size()) {
          return std::unique_ptr<Element>(new ElementAdapter<Scalar<TensorQuad<QuadP1>>>(options));
        }
          /* If elastic fields detected, return a coupled elastic element. */
        else if (el2d) {
          return std::unique_ptr<Element>(new ElementAdapter<ElasticToAcoustic2D<Scalar<TensorQuad<
              QuadP1>>>>(options));
        }
          /* If a mesh boundary is detected but no coupling. */
        else if (bnd && !el2d) {
          return std::unique_ptr<Element>(new ElementAdapter<HomogeneousDirichlet<Scalar<TensorQuad<
              QuadP1>>>>(options));
        }
          /* If a mesh boundary and coupling is detected. */
        else if (bnd && el2d) {
          return std::unique_ptr<Element>(new ElementAdapter<HomogeneousDirichlet<
              ElasticToAcoustic2D<Scalar<TensorQuad<QuadP1>>>>>(options));
        }

      } else if (physics_base == elastic_2d) {
        /* If only elastic, return a base elastic. */
        if (!physics_couple.size()) {
          return std::unique_ptr<Element>(new ElementAdapter<Elastic2D<TensorQuad<QuadP1>>>(options));
        }
          /* If acoustic fields are detected, return a coupled acoustic element. */
        else if (fl) {
          return std::unique_ptr<Element>(new ElementAdapter<AcousticToElastic2D<Elastic2D<
              TensorQuad<QuadP1>>>>(options));
        }
          /* If a mesh boundary is detected, but no coupling. */
        else if (bnd && !fl) {
          return std::unique_ptr<Element>(new ElementAdapter<HomogeneousDirichlet<Elastic2D<
              TensorQuad<QuadP1>>>>(options));
        }
          /* If a mesh boundary and coupling is detected. */
        else if (bnd && fl) {
          return std::unique_ptr<Element>(new ElementAdapter<HomogeneousDirichlet<
              AcousticToElastic2D<Elastic2D<TensorQuad<QuadP1>>>>>(options));
        }
        else {
          throw std::runtime_error(
              "Runtime error: Element physics " + options->PhysicsSystem() + " and element "
                  "type " + options->ElementShape() + " not supported.");
        }
      }
    }
    else if (options->ElementShape() == "hex") {
      /* If only acoustic, return a base acoustic. */
      if (physics_base == fluid) {
        if (!physics_couple.size()) {
          return std::unique_ptr<Element> (
              new ElementAdapter<
                  Scalar<
                      Hexahedra<
                          HexP1>>>(options)); }
        /* If a mesh boundary is detected, but no coupling. */
        else if (bnd && !fl) {
          return std::unique_ptr<Element> (
              new ElementAdapter<
                  HomogeneousDirichlet<
                      Scalar<
                          Hexahedra<
                              HexP1>>>>(options)); }
      } else {
        throw std::runtime_error(
            "Runtime Error: Element physics " + options->PhysicsSystem() + " not supported.");
      }
    }

  } catch (std::exception &e) {
    LOG() << e.what();
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }

  return NULL;

}
