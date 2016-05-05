//
// Created by Michael Afanasiev on 2016-05-03.
//

#include <Element/HyperCube/Quad/QuadP1.h>
#include <Element/HyperCube/QuadNew.h>
#include "AcousticNew.h"

using namespace Eigen;

template <typename Derived>
std::vector<std::string> AcousticNew<Derived>::PullElementalFields() const { return { "u" }; }

template <typename Derived>
std::vector<std::string> AcousticNew<Derived>::PushElementalFields() const { return { "a" }; }

template <typename Derived>
MatrixXd AcousticNew<Derived>::computeStiffnessTerm(
    const Eigen::Ref<const Eigen::MatrixXd> &displacement) {
}

template class AcousticNew<QuadP1>;