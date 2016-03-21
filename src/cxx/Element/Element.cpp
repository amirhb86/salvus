#include "Element.h"

#include <Eigen/Dense>

/*
 * STATIC variables WHICH ARE ONLY ON THE REFERENCE ELEMENT.
 */
int Element2D::mNumberDofVertex;
int Element2D::mNumberDofEdge;
int Element2D::mNumberDofFace;
int Element2D::mNumberIntegrationPoints;
int Element2D::mPolynomialOrder;
Eigen::VectorXi Element2D::mClosureMapping;
Eigen::MatrixXd Element2D::mGradientOperator;
