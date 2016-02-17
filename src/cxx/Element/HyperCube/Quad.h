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
#include "Quad/Autogen/order4_square.h"
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

class Quad {

    // Static functions which only required the reference element.
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

    // Functions which require element data, but do not depend on physical systems.
    bool mCheckHull(double x, double z);
    Eigen::Vector2d inverseCoordinateTransform(const double &x_real, const double &z_real,
                                               double eps, double eta);

protected:

    static const int mNumberVertex = 4;
    static const int mNumberDimensions = 2;

    // Functions which the derived class absolutely needs, but only depends on the reference element.
    static Eigen::Vector4d interpolateShapeFunctions(const double &eps, const double &eta);
    static Eigen::VectorXi mClosureMapping;
    static Eigen::MatrixXd mGradientOperator;
    static Eigen::VectorXd mIntegrationWeightsEps;
    static Eigen::VectorXd mIntegrationWeightsEta;
    static Eigen::VectorXd mIntegrationCoordinatesEps;
    static Eigen::VectorXd mIntegrationCoordinatesEta;
    static int mNumberDofVertex, mNumberDofEdge, mNumberDofFace, mNumberDofVolume;
    static int mNumberIntegrationPointsEps, mNumberIntegrationPointsEta, mNumberIntegrationPoints, mPolynomialOrder;
    static Eigen::Map<const Eigen::VectorXd> epsVectorStride(
            const Eigen::VectorXd &function, const int &eta_index);
    static Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>> etaVectorStride(
            const Eigen::VectorXd &function, const int &eta_index);

    int mElementNumber;
    double mTime;

    std::vector<Source*> mSources;
    Eigen::VectorXd mMassMatrix;
    Eigen::Matrix<double,2,4> mVertexCoordinates;

    Eigen::Matrix<double,2,2> jacobianAtPoint(PetscReal eps, PetscReal eta);

    Eigen::Vector4d __interpolateMaterialProperties(ExodusModel *model,
                                                    std::string parameter_name);


public:

    Quad(Options options);
    virtual Quad * clone() const = 0;

    static Eigen::VectorXd GllPointsForOrder(const int order);
    static Eigen::VectorXd GllIntegrationWeightForOrder(const int order);
    static Eigen::VectorXi ClosureMapping(const int order, const int dimension);

    void scatterMassMatrix(Mesh *mesh);
    void readOperators();
    void attachVertexCoordinates(DM &distributed_mesh);
    void attachSource(std::vector<Source*> sources);

    /**
     * Simple function to set the (remembered) element number.
     */
    void SetLocalElementNumber(const int &element_number) { mElementNumber = element_number; }

    /**
     * Sets the current simulation time. This is used internally, for example, by any sources residing on the element.
     */
    void SetTime(const double &time) { mTime = time; }

    // Attribute gets.
    int NumberDofEdge() const { return mNumberDofEdge; }
    int NumberDofFace() const { return mNumberDofFace; }
    int NumberDofVertex() const { return mNumberDofVertex; }
    int NumberDofVolume() const { return mNumberDofVolume; }
    int NumberDimensions() const { return mNumberDimensions; }

    virtual Eigen::MatrixXd checkOutField(Mesh *mesh, const std::string name);
    virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement) = 0;

    // Pure virtual methods.
    virtual void checkInField(Mesh *mesh) = 0;

    virtual void computeSourceTerm() = 0;
    virtual void computeSurfaceTerm() = 0;

    virtual void assembleMassMatrix() = 0;
    virtual void interpolateMaterialProperties(ExodusModel *model) = 0;

};


#endif //SALVUS_HYPERCUBE_H
