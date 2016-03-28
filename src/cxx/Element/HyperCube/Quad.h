//
// Created by Michael Afanasiev on 2016-01-30.
//

#pragma once

#include <Eigen/Dense>
#include <Model/ExodusModel.h>
#include <Source/Source.h>
#include <Mesh/Mesh.h>
#include <Element/Element.h>

extern "C" {
#include "Quad/Autogen/quad_autogen.h"
};

/**
 * Base class of an abstract four node quadrilateral. The reference element is set up as below.
 *
 * (n3)______________(n2)
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

class Quad : public Element2D {

    /**
     * Shape function contribution from node zero.
     * @param [in] eps Epsilon position in reference quad [-1, 1]
     * @param [in] eta Eta position in reference quad [-1, 1]
     * @returns n0 evaluated at (eps, eta)
     */
    static double n0(const double &eps, const double &eta);

    /**
     * Shape function contribution from node one.
     * @param [in] eps Epsilon position in reference quad [-1, 1]
     * @param [in] eta Eta position in reference quad [-1, 1]
     * @returns n1 evaluated at (eps, eta)
     */
    static double n1(const double &eps, const double &eta);

    /**
     * Shape function contribution from node two.
     * @param [in] eps Epsilon position in reference quad [-1, 1]
     * @param [in] eta Eta position in reference quad [-1, 1]
     * @returns n2 evaluated at (eps, eta)
     */
    static double n2(const double &eps, const double &eta);

    /**
     * Shape function contribution from node three.
     * @param [in] eps Epsilon position in reference quad [-1, 1]
     * @param [in] eta Eta position in reference quad [-1, 1]
     * @returns n3 evaluated at (eps, eta)
     */
    static double n3(const double &eps, const double &eta);

    /**
     * Derivative of node zero's shape function in the eps direction.
     * @param [in] eta Eta position in reference quad [-1, 1]
     * @returns Value at (eps, eta)
     */
    static double dn0deps(const PetscReal &eta);

    /**
     * Derivative of node one's shape function in the eps direction.
     * @param [in] eta Eta position in reference quad [-1, 1]
     * @returns Value at (eps, eta)
     */
    static double dn1deps(const PetscReal &eta);

    /**
     * Derivative of node two's shape function in the eps direction.
     * @param [in] eta Eta position in reference quad [-1, 1]
     * @returns Value at (eps, eta)
     */
    static double dn2deps(const PetscReal &eta);

    /**
     * Derivative of node three's shape function in the eps direction.
     * @param [in] eta Eta position in reference quad [-1, 1]
     * @returns Value at (eps, eta)
     */
    static double dn3deps(const PetscReal &eta);

    /**
     * Derivative of node zero's shape function in the eta direction.
     * @param [in] eta Eta position in reference quad [-1, 1]
     * @returns Value at (eps, eta)
     */
    static double dn0deta(const PetscReal &eps);

    /**
     * Derivative of node one's shape function in the eta direction.
     * @param [in] eta Eta position in reference quad [-1, 1]
     * @returns Value at (eps, eta)
     */
    static double dn1deta(const PetscReal &eps);

    /**
     * Derivative of node two's shape function in the eta direction.
     * @param [in] eta Eta position in reference quad [-1, 1]
     * @returns Value at (eps, eta)
     */
    static double dn2deta(const PetscReal &eps);

    /**
     * Derivative of node three's shape function in the eta direction.
     * @param [in] eta Eta position in reference quad [-1, 1]
     * @returns Value at (eps, eta)
     */
    static double dn3deta(const PetscReal &eps);
    
protected:

    /*****************************************************************************
     * STATIC MEMBERS. THESE VARIABLES AND FUNCTIONS SHOULD APPLY TO ALL ELEMENTS.
     *****************************************************************************/

    static int mNumberVertex;         /** < Number of element vertices. */
    
    static int mNumberIntegrationPointsEps;     /** < Number of integration points in the epsilon direction (e.g. 5 for a 4th order gll basis) */
    static int mNumberIntegrationPointsEta;     /** < Number of integration points in the eta direction (e.g. 5 for a 4th order gll basis) */

    static Eigen::MatrixXd mGradientOperator;           /** < Derivative of shape function n (col) at pos. m (row) */
    
    static Eigen::VectorXd mIntegrationWeightsEps;      /** < Integration weights along epsilon direction. */
    static Eigen::VectorXd mIntegrationWeightsEta;      /** < Integration weights along eta direction. */
    static Eigen::VectorXd mIntegrationCoordinatesEps;  /** < Integration points along epsilon direction */
    static Eigen::VectorXd mIntegrationCoordinatesEta;  /** < Integration points along eta direction */

    /**
     * Returns the shape function coefficients for a given location (eps, eta) in the reference cube.
     * @param [in] eps Epsilon in the reference element.
     * @param [in] eta Eta in the reference element.
     * @returns A vector containing the [4] coefficients from each shape function.
     */
    static Eigen::Vector4d interpolateShapeFunctions(const double &eps, const double &eta);

    /**
     * Returns the lagrange polynomial coefficients for a given location (eps, eta) in the reference cube.
     * @param [in] eps Epsilon on the reference element.
     * @param [in] eta Eta on the reference element.
     * @returns A vector containing the polynomial coefficient at each gll point.
     */
    static Eigen::VectorXd interpolateLagrangePolynomials(const double eps, const double eta, const int p_order);

    /**
     * Sets up the proper stride for a field in the epsilon direction, on the tensorized gll basis.
     * In many instances, the only non-zero contributions to integrals or derivatives on the gll basis are those
     * contributions which lie along a specific direction in the reference element. Given a properly ordered field
     * defined on all the gll points, this function returns a `const` view into that vector. This view includes all
     * points of the given field along a certain 'row' (as we are working in the epsilon here), from the beginning to
     * end of that 'row'.
     *
     * As an example, consider the 4th order gll basis. If a field u is defined in an ordered manner across all 25
     * gll points, a call to epsVectorStride(u, 1) would return the field at indices 5, 6, 7, 8, and 9 -- the row
     * which the eta_index refers to.
     *
     * @param [in] function Field defined across an element, for which a view is desired.
     * @param [in] eta_index Eta_index in the reference element -- to determine which epsilon 'row' we're on.
     * @returns A const pointer with the proper stride and starting point for the desired field points.
     */
    static Eigen::Map<const Eigen::VectorXd> epsVectorStride(
            const Eigen::VectorXd &function, const int &eta_index);

    /**
     * Sets up the proper stride for a field in the eta direction, on the tensorized gll basis.
     * In many instances, the only non-zero contributions to integrals or derivatives on the gll basis are those
     * contributions which lie along a specific direction in the reference element. Given a properly ordered field
     * defined on all the gll points, this function returns a `const` view into that vector. This view includes all
     * points of the given field along a certain 'column' (as we are working in the eta here), from the beginning to
     * end of that 'column'.
     *
     * As an example, consider the 4th order gll basis. If a field u is defined in an ordered manner across all 25
     * gll points, a call to etaVectorStride(u, 1) would return the field at indices 1, 6, 11, 16, and 21 -- the column
     * which the eta_index refers to.
     *
     * @param [in] function Field defined across an element, for which a view is desired.
     * @param [in] eps_index Eps_index in the reference element -- to determine which eta 'column' we're on.
     * @returns A const pointer with the proper stride and starting point for the desired field points.
     */
    static Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<>> etaVectorStride(
            const Eigen::VectorXd &function, const int &eps_index);

    /**********************************************************************************
     * OBJECT MEMBERS. THESE VARIABLES AND FUNCTIONS SHOULD APPLY TO SPECIFIC ELEMENTS.
    ***********************************************************************************/

public:
    
    /**
     * Constructor.
     * Sets quantities such as number of dofs, among other things, from the options class.
     * @param [in] options Populated options class.
     */
    Quad(Options options);

    
    /**
     * Returns the quadrature locations for a given polynomial order.
     * @param [in] order The polynmomial order.
     * @returns Vector of GLL points.
     */
    static Eigen::VectorXd GllPointsForOrder(const int order);

    /**
     * Returns the quadrature intergration weights for a polynomial order.
     * @param [in] order The polynomial order.
     * @returns Vector of quadrature weights.
     */
    static Eigen::VectorXd GllIntegrationWeightForOrder(const int order);

    /**
     * Returns the mapping from the PETSc to Salvus closure.
     * @param [in] order The polynomial order.
     * @param [in] dimension Element dimension.
     * @returns Vector containing the closure mapping (field(closure(i)) = petscField(i))
     */
    static Eigen::VectorXi ClosureMapping(const int order, const int dimension);

    // face closure not currently possible...
    // static Eigen::VectorXi GetFaceClosureMapping() { return mFaceClosureMapping; }

    /********************************************************
     ********* Inherited methods from Element 2D ************
     ********************************************************/

    /**
     * Checks whether a given point in realspace (x, z) is within the current element.
     * A simple convex hull algorithm is implemented. Should work as long as the sides of the element are straight
     * lines, but will likely fail for higher order shape functions.
     * @param [in] x X-coordinate in real space.
     * @param [in] z Z-coordinate in real space.
     * @returns True if point is inside, False if not.
     */
    bool mCheckHull(double x, double z);

    /**
     * Given a point in realspace, determines the equivalent location in the reference element.
     * Since the shape function are bilinear, the cross terms introduce nonlinearities into the shape function
     * interpolation. As such, we used a simple implementation of Newton's method to minimize the objective function:
     * z = x_real - shape_functions * vertex_coordinates(eta, eps), for each coordinate, where (eps, eta) are the
     * primal variables.
     * @param [in] x_real X-coordinate in real space.
     * @param [in] z_real Z-coordinate in real space.
     * @return A Vector (eps, eta) containing the coordinates in the reference element.
     */
    Eigen::Vector2d inverseCoordinateTransform(const double &x_real, const double &z_real);
    
    /**
     * 2x2 Jacobian matrix at a point (eps, eta).
     * This method returns an Eigen::Matrix representation of the Jacobian at a particular point. Since the return
     * value is this Eigen object, we can perform additional tasks on the return value (such as invert, and get
     * determinant).
     * @param [in] eps Epsilon position on the reference element.
     * @param [in] eta Eta position on the reference element.
     * @returns 2x2 statically allocated Eigen matrix object.
     */
    Eigen::Matrix<double,2,2> jacobianAtPoint(PetscReal eps, PetscReal eta);

    /**
     * 2x2 inverse Jacobian matrix at a point (eps, eta).  This method returns an Eigen::Matrix
     * representation of the inverse Jacobian at a particular point.
     * @param [in] eps Epsilon position on the reference element.
     * @param [in] eta Eta position on the reference element.     
     * @returns (inverse Jacobian matrix,determinant of that matrix) as a `std::tuple`. Tuples can be
     * "destructured" using a `std::tie`.
     */
    std::tuple<Eigen::Matrix2d,PetscReal> inverseJacobianAtPoint(PetscReal eps, PetscReal eta);

    /**
     * Attaches a material parameter to the vertices on the current element.
     * Given an exodus model object, we use a kD-tree to find the closest parameter to a vertex. In practice, this
     * closest parameter is often exactly coincident with the vertex, as we use the same exodus model for our mesh
     * as we do for our parameters.
     * @param [in] model An exodus model object.
     * @param [in] parameter_name The name of the field to be added (i.e. velocity, c11).
     * @returns A Vector with 4-entries... one for each Element vertex, in the ordering described above.
     */
    Eigen::Vector4d __interpolateMaterialProperties(ExodusModel *model,
                                                            std::string parameter_name);

    /**
     * Utility function to integrate a field over the element. This could probably be made static, but for now I'm
     * just using it to check some values.
     * @param [in] field The field which to integrate, defined on each of the gll points.
     * @returns The scalar value of the field integrated over the element.
     */
    double integrateField(const Eigen::VectorXd &field);

    /**
     * Setup the auto-generated gradient operator, and stores the result in mGradientOperator.
     */
    void setupGradientOperator();

    
    /**
     * Queries the passed DM for the vertex coordinates of the specific element. These coordinates are saved
     * in mVertexCoordiantes.
     * @param [in] distributed_mesh PETSc DM object.
     *
     */
    void attachVertexCoordinates(DM &distributed_mesh);

    /**
     * Attach source.
     * Given a vector of abstract source objects, this function will query each for its spatial location. After
     * performing a convex hull test, it will perform a quick inverse problem to determine the position of any sources
     * within each element in reference coordinates. These reference coordinates are then saved in the source object.
     * References to any sources which lie within the element are saved in the mSources vector.
     * @param [in] sources A vector of all the sources defined for a simulation run.
     */
    void attachSource(std::vector<Source*> sources);

    /**
     * Simple function to set the (remembered) element number.
     */
    void SetLocalElementNumber(int element_number) { mElementNumber = element_number; }
           
    /**
     * Builds nodal coordinates (x,z) on all mesh degrees of freedom.
     * @param mesh [in] The mesh.
     */
    std::tuple<Eigen::VectorXd,Eigen::VectorXd> buildNodalPoints(Mesh* mesh);

};
