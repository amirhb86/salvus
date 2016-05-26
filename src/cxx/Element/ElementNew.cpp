
#include <ElementNew.h>
#include <ElementAdapter.h>

#include <HyperCube/QuadNew.h>
#include <HyperCube/Quad/QuadP1.h>

#include <Simplex/TriangleNew.h>
#include <Simplex/Triangle/TriP1.h>

#include <Simplex/TetrahedraNew.h>
#include <Simplex/Tetrahedra/TetP1.h>

#include <HyperCube/HexahedraNew.h>
#include <HyperCube/Hex/HexP1.h>

#include <Physics/AcousticNew.h>
#include <Physics/AcousticTriNew.h>
#include <Physics/Acoustic3D.h>
#include <Physics/AcousticTet3D.h>
#include <Physics/Acoustic3D_V.h>
#include <Physics/ElasticNew.h>

/* Define all possible element classes as types here. */
typedef class ElementAdapter<AcousticNew<QuadNew<QuadP1>>> AcousticQuadP1;
typedef class ElementAdapter<AcousticTriNew<TriangleNew<TriP1>>> AcousticTriP1;
typedef class ElementAdapter<Acoustic3D<HexahedraNew<HexP1>>> AcousticHexP1;
typedef class ElementAdapter<Acoustic3D<TetrahedraNew<TetP1>>> AcousticTetP1v2;
typedef class ElementAdapter<Acoustic3D_V<HexahedraNew<HexP1>>> AcousticVHexP1;
typedef class ElementAdapter<AcousticTet3D<TetrahedraNew<TetP1>>> AcousticTetP1;
typedef class ElementAdapter<ElasticNew<QuadNew<QuadP1>>> ElasticQuadP1;

std::shared_ptr<ElementNew> ElementNew::Factory(Options options) {
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
