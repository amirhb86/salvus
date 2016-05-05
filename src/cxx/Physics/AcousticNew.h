#pragma once

#include <Eigen/Dense>
#include <Element/Element.h>
#include <vector>

template <typename Shape>
class AcousticNew: public Shape {

 public:

  // Delegates.


//  std::shared_ptr<Element> clone() const { return std::shared_ptr<Element> (new AcousticNew(*this)); }
  std::vector<std::string> PullElementalFields() const;
  std::vector<std::string> PushElementalFields() const;

  Eigen::MatrixXd computeStiffnessTerm(const Eigen::Ref<const Eigen::MatrixXd>& displacement);

};

