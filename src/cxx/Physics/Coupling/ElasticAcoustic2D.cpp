#include <Physics/ElasticAcoustic2D.h>
#include <Utilities/Options.h>

template <typename BasePhysics>
ElasticAcoustic2D<BasePhysics>::ElasticAcoustic2D(Options options): BasePhysics(options) {

}



#include <Physics/Elastic2D.h>
#include <Element/HyperCube/Quad.h>
#include <Element/HyperCube/QuadP1.h>
template class ElasticAcoustic2D<Elastic2D<Quad<QuadP1>>>;
