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
#include <Physics/AcousticTri.h>
#include <Physics/Acoustic3D.h>
#include <Physics/AcousticHex3D.h>
#include <Physics/AcousticHex3D_LF.h>
#include <Physics/AcousticTet.h>
#include <Physics/Acoustic3D_V.h>
#include <Physics/Elastic2D.h>
#include <Physics/Elastic3D.h>
#include <Physics/AcousticElastic2D.h>
#include <Physics/ElasticAcoustic2D.h>
#include <Utilities/Utilities.h>

/* Define all possible element classes as types here. */
typedef class ElementAdapter<Scalar<TensorQuad<QuadP1>>> AcousticQuadP1;
typedef class ElementAdapter<AcousticTri<Triangle<TriP1>>> AcousticTriP1;
typedef class ElementAdapter<Acoustic3D<Hexahedra<HexP1>>> AcousticHexP1;
typedef class ElementAdapter<AcousticHex3D<Hexahedra<HexP1>>> AcousticHexP1v2;
typedef class ElementAdapter<AcousticHex3D_LF<Hexahedra<HexP1>>> AcousticHexP1_fast_lf;
typedef class ElementAdapter<Acoustic3D<Tetrahedra<TetP1>>> AcousticTetP1v2;
typedef class ElementAdapter<Acoustic3D_V<Hexahedra<HexP1>>> AcousticVHexP1;
typedef class ElementAdapter<AcousticTet<Tetrahedra<TetP1>>> AcousticTetP1;
typedef class ElementAdapter<Elastic2D<TensorQuad<QuadP1>>> ElasticQuadP1;
typedef class ElementAdapter<Elastic3D<Hexahedra<HexP1>>> ElasticHexP1;

/* Coupled classes. */
typedef class ElementAdapter<AcousticToElastic2D<Elastic2D<TensorQuad<QuadP1>>>> AcousticCplElasticQuadP1;
typedef class ElementAdapter<ElasticToAcoustic2D<Scalar<TensorQuad<QuadP1>>>> ElasticCplAcousticQuadP1;

std::unique_ptr<Element> Element::Factory(const std::vector<std::string>& physics_base,
                                          const std::vector<std::string>& physics_couple,
                                          std::unique_ptr<Options> const &options) {

  // define field combinations.
  std::vector<std::string> acoustic_fields = {"u"};
  std::vector<std::string> elastic_2d_fields = {"ux", "uy"};
  std::vector<std::string> elastic_3d_fields = {"ux", "uy", "uz"};

  try {
    if (options->ElementShape() == "quad_new") {
      if (physics_base == acoustic_fields) {

        /* If only acoustic, return a base acoustic. */
        if (!physics_couple.size()) {
          return std::unique_ptr<Element> (new AcousticQuadP1(options));
        }

        /* If elastic fields detected, return a coupled elastic element. */
        else if (physics_couple == elastic_2d_fields) {
          return std::unique_ptr<Element> (new ElasticCplAcousticQuadP1(options));
        }

      } else if (physics_base == elastic_2d_fields) {

        /* If only elastic, return a base elastic. */
        if (!physics_couple.size()) {
          return std::unique_ptr<Element>(new ElasticQuadP1(options));
        }

        /* If acoustic fields are detected, return a coupled acoustic element. */
        else if (physics_couple == acoustic_fields) {
          return std::unique_ptr<Element>(new AcousticCplElasticQuadP1(options));
        }

      } else {
        throw std::runtime_error("Runtime Error: Element physics " + options->PhysicsSystem() + " not supported.");
      }
    }

    else if (options->ElementShape() == "triangle_new") {
      if (options->PhysicsSystem() == "acoustic") {
        return std::unique_ptr<Element>(new AcousticTriP1(options));
      } else {
        throw std::runtime_error("Runtime Error: Element physics " + options->PhysicsSystem() + " not supported.");
      }
    }
    else if (options->ElementShape() == "hex_new") {
      if (physics_base == acoustic_fields) {
        if (!physics_couple.size()) {
          return std::unique_ptr<Element>(new AcousticVHexP1(options));
        }
      } else if (physics_base == elastic_3d_fields) {
        if (!physics_couple.size()) {
          return std::unique_ptr<Element>(new ElasticHexP1(options));
        }
      }

      if (options->PhysicsSystem() == "acoustic") {
        return std::unique_ptr<Element> (new AcousticHexP1(options));
      } else if (options->PhysicsSystem() == "acoustic_fast") {
        return std::unique_ptr<Element> (new AcousticHexP1v2(options));
      } else if (options->PhysicsSystem() == "acoustic_v") {
        std::cout << "HI" << std::endl;
        return std::unique_ptr<Element> (new AcousticVHexP1(options));
      } else if (options->PhysicsSystem() == "acoustic_lf") {
        return std::unique_ptr<Element> (new AcousticHexP1_fast_lf(options));
      } else {
        throw std::runtime_error("Runtime Error: Element physics " + options->PhysicsSystem() + " not supported.");
      }
    } else if (options->ElementShape() == "tet_new") {
      if (options->PhysicsSystem() == "acoustic") {
        return std::unique_ptr<Element>(new AcousticTetP1(options));
      } else {
        throw std::runtime_error("Runtime Error: Element physics " + options->PhysicsSystem() + " not supported.");
      }
    } else {
      throw std::runtime_error("Runtime Error: Element shape " + options->ElementShape() + " not supported.");
    }
  } catch (std::exception &e) {
    std::cout << e.what();
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }

  return NULL;
  MPI_Abort(PETSC_COMM_WORLD, -1);
}
