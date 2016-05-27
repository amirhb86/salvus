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
#include <Physics/AcousticTet.h>
#include <Physics/Acoustic3D_V.h>
#include <Physics/Elastic2D.h>
#include <Utilities/Utilities.h>

/* Define all possible element classes as types here. */
typedef class ElementAdapter<Acoustic2D<Quad<QuadP1>>> AcousticQuadP1;
typedef class ElementAdapter<AcousticTri<Triangle<TriP1>>> AcousticTriP1;
typedef class ElementAdapter<Acoustic3D<Hexahedra<HexP1>>> AcousticHexP1;
typedef class ElementAdapter<Acoustic3D<Tetrahedra<TetP1>>> AcousticTetP1v2;
typedef class ElementAdapter<Acoustic3D_V<Hexahedra<HexP1>>> AcousticVHexP1;
typedef class ElementAdapter<AcousticTet<Tetrahedra<TetP1>>> AcousticTetP1;
typedef class ElementAdapter<Elastic2D<Quad<QuadP1>>> ElasticQuadP1;

std::shared_ptr<Element> Element::Factory(Options options) {
  try {
    if (options.ElementShape() == "quad_new") {
      if (options.PhysicsSystem() == "acoustic") {
        return std::make_shared<AcousticQuadP1>(options);
      } else if (options.PhysicsSystem() == "elastic") {
        return std::make_shared<ElasticQuadP1>(options);
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
      } else if (options.PhysicsSystem() == "acoustic_v") {
        return std::make_shared<AcousticVHexP1>(options);
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

  return nullptr;
}
