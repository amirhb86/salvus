#include <Mesh/Mesh.h>
#include <Source/Source.h>
#include <Physics/Scalar.h>
#include <Physics/ScalarTri.h>
#include <Utilities/Options.h>
#include <Utilities/Logging.h>
#include <Model/ExodusModel.h>

using namespace Eigen;

template <typename Element>
RealMat ScalarTri<Element>::computeStiffnessTerm(const Ref<const RealMat>& u) {
  
  mStiff = Element::StiffnessMatrix()*u.col(0);
  return mStiff;

}

#include <Element/Simplex/Triangle.h>
#include <Element/Simplex/TriP1.h>

template class ScalarTri<Scalar<Triangle<TriP1>>>;

