#include <Physics/AcousticElastic2D.h>
#include <Utilities/Options.h>

template <typename BasePhysics>
AcousticElastic2D<BasePhysics>::AcousticElastic2D(Options options): BasePhysics(options) {

}



#include <Physics/Acoustic2D.h>
#include <Element/HyperCube/Quad.h>
#include <Element/HyperCube/QuadP1.h>
template class AcousticElastic2D<Acoustic2D<Quad<QuadP1>>>;