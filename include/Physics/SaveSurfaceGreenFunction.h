#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <Utilities/Types.h>

template <typename Element>
class SaveSurfaceGreenFunction: public Element {

 public:

  SaveSurfaceGreenFunction<Element>(std::unique_ptr<Options> const &options);

};

