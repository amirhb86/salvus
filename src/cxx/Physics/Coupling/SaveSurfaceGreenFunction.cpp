#include <Physics/SaveSurfaceGreenFunction.h>

template <typename Element>
SaveSurfaceGreenFunction<Element>::SaveSurfaceGreenFunction(
    std::unique_ptr<Options> const &options) : Element(options) {

  std::cout << "HI!" << std::endl;

}

#include <Physics/Scalar.h>
#include <Element/HyperCube/TensorQuad.h>
#include <Element/HyperCube/QuadP1.h>

template class SaveSurfaceGreenFunction<
    Scalar<
        TensorQuad<
            QuadP1>>>;
