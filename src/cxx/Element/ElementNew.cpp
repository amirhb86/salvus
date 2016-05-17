
#include <ElementNew.h>
#include <ElementAdapter.h>

#include <HyperCube/QuadNew.h>
#include <HyperCube/Quad/QuadP1.h>

#include <Physics/AcousticNew.h>
#include <Physics/ElasticNew.h>

/* Define all possible element classes as types here. */
typedef class ElementAdapter<AcousticNew<QuadNew<QuadP1>>> AcousticQuadP1;
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
    } else {
      throw std::runtime_error("Runtime Error: Element shape " + options.ElementShape() + " not supported.");
    }
  } catch (std::exception &e) {
    PRINT_ROOT() << e.what();
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }

  return nullptr;
}