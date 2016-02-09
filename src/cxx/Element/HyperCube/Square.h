//
// Created by Michael Afanasiev on 2016-01-30.
//

#ifndef SALVUS_HYPERCUBE_H
#define SALVUS_HYPERCUBE_H

#include <Eigen/Dense>
#include "../../Model/ExodusModel.h"
#include "../../Source.h"
#include "../../Mesh.h"

extern "C" {
#include "Square/Autogen/order4_square.h"
};

/*
 * Base class of an abstract four node square. The reference element is set up as below.
 *
 * (n2)______________(n3)
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * |                    |
 * (n0)______________(n1)
 *
 * (eta)
 *   ^
 *   |
 *   |
 *   |______> (eps)
*/

class Square {

    static int element_increment;

    static double n0(const double &eps, const double &eta);
    static double n1(const double &eps, const double &eta);
    static double n2(const double &eps, const double &eta);
    static double n3(const double &eps, const double &eta);

    static double dn0deps(const PetscReal &eta);
    static double dn1deps(const PetscReal &eta);
    static double dn2deps(const PetscReal &eta);
    static double dn3deps(const PetscReal &eta);

    static double dn0deta(const PetscReal &eps);
    static double dn1deta(const PetscReal &eps);
    static double dn2deta(const PetscReal &eps);
    static double dn3deta(const PetscReal &eps);

    bool mCheckHull(double x, double z);
    Eigen::Vector2d inverseCoordinateTransform(const double &x_real, const double &z_real,
                                               double eps, double eta);

protected:

    static const int mNumberVertex = 4;
    static const int mNumberDimensions = 2;

    bool mContainsSource;

    int mElementNumber;
    int mNumberDofVertex, mNumberDofEdge, mNumberDofFace, mNumberDofVolume;
    int mNumberIntegrationPointsEps, mNumberIntegrationPointsEta, mNumberIntegrationPoints, mPolynomialOrder;

    double mTime;

    std::vector<Source*> mSources;
    Eigen::VectorXi mClosureMapping;
    Eigen::MatrixXd mGradientOperator;
    Eigen::VectorXd mIntegrationWeightsEps;
    Eigen::VectorXd mIntegrationWeightsEta;
    Eigen::VectorXd mIntegrationCoordinatesEps;
    Eigen::VectorXd mIntegrationCoordinatesEta;
    Eigen::Matrix<double,2,4> mVertexCoordinates;

    Eigen::Matrix<double,2,2> jacobianAtPoint(PetscReal eps, PetscReal eta);
    Eigen::Vector4d interpolateShapeFunctions(PetscReal eps, PetscReal eta);

    Eigen::Vector4d __interpolateMaterialProperties(ExodusModel &model,
                                                    std::string parameter_name);

    Eigen::Map<Eigen::VectorXd> epsVectorStride(
            Eigen::VectorXd &function, const int &eta_index);
    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<>> etaVectorStride(
            Eigen::VectorXd &function, const int &eta_index);


public:

    Square(Options options);
    virtual Square* clone() const = 0;

    static Eigen::VectorXd GllPointsForOrder(const int order);
    static Eigen::VectorXd GllIntegrationWeightForOrder(const int order);
    static Eigen::VectorXi ClosureMapping(const int order, const int dimension);

    void readOperators();
    void attachVertexCoordinates(DM &distributed_mesh);
    virtual void attachSource(std::vector<Source*> sources);

    // Attribute sets.
    void SetLocalElementNumber(const int &element_number) { mElementNumber = element_number; }

    // Attribute gets.
    int NumberDofEdge() const { return mNumberDofEdge; }
    int NumberDofFace() const { return mNumberDofFace; }
    int NumberDofVertex() const { return mNumberDofVertex; }
    int NumberDofVolume() const { return mNumberDofVolume; }
    int NumberDimensions() const { return mNumberDimensions; }

    // Pure virtual methods.
    virtual void registerFieldVectors(Mesh *mesh) = 0;
    virtual void constructStiffnessMatrix(Mesh *mesh) = 0;
    virtual void interpolateMaterialProperties(ExodusModel &model) = 0;

};


#endif //SALVUS_HYPERCUBE_H
