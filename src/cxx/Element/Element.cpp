#include <stdexcept>

#include <Element/Element.h>
#include <Element/ElementAdapter.h>

#include <Element/HyperCube/Quad.h>
#include <Element/HyperCube/QuadP1.h>
#include <Element/Simplex/Triangle.h>
#include <Element/Simplex/TriP1.h>
#include <Element/Simplex/Tetrahedra.h>
#include <Element/Simplex/TetP1.h>
#include <Element/HyperCube/Hexahedra.h>
#include <Element/HyperCube/HexP1.h>

#include <Physics/Acoustic2D.h>
#include <Physics/AcousticTri.h>
#include <Physics/Acoustic3D.h>
#include <Physics/AcousticHex3D.h>
#include <Physics/AcousticHex3D_LF.h>
#include <Physics/AcousticTet.h>
#include <Physics/Acoustic3D_V.h>
#include <Physics/Elastic2D.h>
#include <Physics/AcousticElastic2D.h>
#include <Physics/ElasticAcoustic2D.h>
#include <Utilities/Utilities.h>
#include <stdexcept>

/* Define all possible element classes as types here. */
typedef class ElementAdapter<Acoustic2D<Quad<QuadP1>>> AcousticQuadP1;
typedef class ElementAdapter<AcousticTri<Triangle<TriP1>>> AcousticTriP1;
typedef class ElementAdapter<Acoustic3D<Hexahedra<HexP1>>> AcousticHexP1;
typedef class ElementAdapter<AcousticHex3D<Hexahedra<HexP1>>> AcousticHexP1v2;
typedef class ElementAdapter<AcousticHex3D_LF<Hexahedra<HexP1>>> AcousticHexP1_fast_lf;
typedef class ElementAdapter<Acoustic3D<Tetrahedra<TetP1>>> AcousticTetP1v2;
typedef class ElementAdapter<Acoustic3D_V<Hexahedra<HexP1>>> AcousticVHexP1;
typedef class ElementAdapter<AcousticTet<Tetrahedra<TetP1>>> AcousticTetP1;
typedef class ElementAdapter<Elastic2D<Quad<QuadP1>>> ElasticQuadP1;

/* Coupled classes. */
typedef class ElementAdapter<AcousticToElastic2D<Elastic2D<Quad<QuadP1>>>> AcousticCplElasticQuadP1;
typedef class ElementAdapter<ElasticToAcoustic2D<Acoustic2D<Quad<QuadP1>>>> ElasticCplAcousticQuadP1;

std::shared_ptr<Element> Element::Factory(const std::vector<std::string>& physics_base,
                                          const std::vector<std::string>& physics_couple,
                                          Options options) {

  // define field combinations.
  std::vector<std::string> acoustic_fields = {"u"};
  std::vector<std::string> elastic_2d_fields = {"ux", "uy"};

  try {
    if (options.ElementShape() == "quad_new") {
      if (physics_base == acoustic_fields) {
        /* If only acoustic, return a base acoustic. */
        if (!physics_couple.size()) {
          return std::make_shared<AcousticQuadP1>(options);
        }
        /* If elastic fields detected, return a coupled elastic element. */
        else if (physics_couple == elastic_2d_fields) {
          return std::make_shared<ElasticCplAcousticQuadP1>(options);
        }
      } else if (physics_base == elastic_2d_fields) {
        /* If only elastic, return a base elastic. */
        if (!physics_couple.size()) {
          return std::make_shared<ElasticQuadP1>(options);
        }
        /* If acoustic fields are detected, return a coupled acoustic element. */
        else if (physics_couple == acoustic_fields) {
          return std::make_shared<AcousticCplElasticQuadP1>(options);
        }
      } else {
        throw std::runtime_error("Runtime Error: Element physics " + options.PhysicsSystem() + " not supported.");
      }
    }
    else if (options.ElementShape() == "triangle_new") {
      if (options.PhysicsSystem() == "acoustic") {
        return std::make_shared<AcousticTriP1>(options);
      } else {
        throw std::runtime_error("Runtime Error: Element physics " + options.PhysicsSystem() + " not supported.");
      }
    }      
    else if (options.ElementShape() == "hex_new") {
      if (options.PhysicsSystem() == "acoustic") {
        return std::make_shared<AcousticHexP1>(options);
      } else if (options.PhysicsSystem() == "acoustic_fast") {
        return std::make_shared<AcousticHexP1v2>(options);
      } else if (options.PhysicsSystem() == "acoustic_v") {
        return std::make_shared<AcousticVHexP1>(options);
      } else if (options.PhysicsSystem() == "acoustic_lf") {
        return std::make_shared<AcousticHexP1_fast_lf>(options);
      } else {
        throw std::runtime_error("Runtime Error: Element physics " + options.PhysicsSystem() + " not supported.");
      }
    } else if (options.ElementShape() == "tet_new") {
      if (options.PhysicsSystem() == "acoustic") {
        return std::make_shared<AcousticTetP1>(options);
      } else {
        throw std::runtime_error("Runtime Error: Element physics " + options.PhysicsSystem() + " not supported.");
      }
    } else {
      throw std::runtime_error("Runtime Error: Element shape " + options.ElementShape() + " not supported.");
    }
  } catch (std::exception &e) {
    PRINT_ROOT() << e.what();
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }

  MPI_Abort(PETSC_COMM_WORLD, -1);
}
