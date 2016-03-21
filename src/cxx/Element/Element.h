#pragma once

#include <Eigen/Dense>
#include <Model/ExodusModel.h>
#include <Source/Source.h>
#include <Mesh/Mesh.h>

/**
 * Base class for 2D elements
 */
class Element2D {
    
protected:
    
    /*****************************************************************************
     * STATIC MEMBERS. THESE VARIABLES AND FUNCTIONS SHOULD APPLY TO ALL ELEMENTS.
     *****************************************************************************/
    static const int mNumberDimensions = 2;     /** < Number of element dimensions. */
    
    static int mNumberIntegrationPoints;        /** < Total number of integration points (e.g. 25 for a 4th order gll basis) */

    static int mNumberDofVertex;                /** < Number of dofs on a vertex (e.g. 1 for 4th order gll basis) */
    static int mNumberDofEdge;                  /** < Number of dofs on an edge (e.g. 3 for 4th order gll basis) */
    static int mNumberDofFace;                  /** < Number of dofs on a face (eg. 9 for a 4th order gll basis) */

    static int mPolynomialOrder;                /** < base Lagrange polynomial order */

    static Eigen::VectorXi mClosureMapping;             /** < Mapping from our element closure
                                                            numbering to PETSc's */

    // not currently possible...
    // static Eigen::VectorXi mFaceClosureMapping;             /** < Mapping from our face closure
    //                                                             numbering to PETSc's */
    
    static Eigen::MatrixXd mGradientOperator;           /** < Derivative of shape function n (col) at pos. m (row) */

    /**********************************************************************************
     * OBJECT MEMBERS. THESE VARIABLES AND FUNCTIONS SHOULD APPLY TO SPECIFIC ELEMENTS.
     ***********************************************************************************/

    int mElementNumber; /** Element number on the local processor. */

    std::vector<Source*> mSources;  /** Vector of abstract sources belonging (spatial) to the element */

    Eigen::VectorXd mMassMatrix;    /** Elemental mass matrix */
    
    Eigen::MatrixXd mVertexCoordinates;   /** Vertex coordinates ordered as above. row(0)->x, row(1)->z */

    Eigen::Matrix<double,2,1> mElementCenter; /** (x, z) location of element center */

    bool mOnBoundary;    /** < Whether or not the current element has a special boundary condition */

    std::map<std::string,int> mBoundaries;  /** < Map relating the type of boundary to the Petsc edge number */

public:

    /**
     * Copy constructor.
     * Returns a copy. Use this once the reference element is set up via the constructor, to allocate space for
     * all the unique elements on a processor.
     */
    virtual Element2D * clone() const = 0;
    
    /**
     * Checks whether a given point in realspace (x, z) is within the current element.
     * A simple convex hull algorithm is implemented. Should work as long as the sides of the element are straight
     * lines, but will likely fail for higher order shape functions.
     * @param [in] x X-coordinate in real space.
     * @param [in] z Z-coordinate in real space.
     * @returns True if point is inside, False if not.
     */
    virtual bool mCheckHull(double x, double z) = 0;
    
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
    virtual Eigen::Vector2d inverseCoordinateTransform(const double &x_real, const double &z_real) = 0;
    
    /**
     * 2x2 Jacobian matrix at a point (eps, eta).
     * This method returns an Eigen::Matrix representation of the Jacobian at a particular point. Since the return
     * value is this Eigen object, we can perform additional tasks on the return value (such as invert, and get
     * determinant).
     * @param [in] eps Epsilon position on the reference element.
     * @param [in] eta Eta position on the reference element.
     * @returns 2x2 statically allocated Eigen matrix object.
     */
    virtual Eigen::Matrix<double,2,2> jacobianAtPoint(PetscReal eps, PetscReal eta) = 0;

    /**
     * 2x2 inverse Jacobian matrix at a point (eps, eta).  This method returns an Eigen::Matrix
     * representation of the inverse Jacobian at a particular point.
     * @param [in] eps Epsilon position on the reference element.
     * @param [in] eta Eta position on the reference element.     
     * @returns (inverse Jacobian matrix,determinant of that matrix) as a `std::tuple`. Tuples can be
     * "destructured" using a `std::tie`.
     */
    virtual std::tuple<Eigen::Matrix2d,PetscReal> inverseJacobianAtPoint(PetscReal eps, PetscReal eta) = 0;

    /**
     * Attaches a material parameter to the vertices on the current element.
     * Given an exodus model object, we use a kD-tree to find the closest parameter to a vertex. In practice, this
     * closest parameter is often exactly coincident with the vertex, as we use the same exodus model for our mesh
     * as we do for our parameters.
     * @param [in] model An exodus model object.
     * @param [in] parameter_name The name of the field to be added (i.e. velocity, c11).
     * @returns A Vector with 4-entries... one for each Element vertex, in the ordering described above.
     */
    virtual Eigen::Vector4d __interpolateMaterialProperties(ExodusModel *model,
                                                            std::string parameter_name) = 0;

    /**
     * Utility function to integrate a field over the element. This could probably be made static, but for now I'm
     * just using it to check some values.
     * @param [in] field The field which to integrate, defined on each of the gll points.
     * @returns The scalar value of the field integrated over the element.
     */
    virtual double integrateField(const Eigen::VectorXd &field) = 0;

    /**
     * Setup the auto-generated gradient operator, and stores the result in mGradientOperator.
     */
    virtual void setupGradientOperator() = 0;

    
    /**
     * Queries the passed DM for the vertex coordinates of the specific element. These coordinates are saved
     * in mVertexCoordiantes.
     * @param [in] distributed_mesh PETSc DM object.
     *
     */
    virtual void attachVertexCoordinates(DM &distributed_mesh) = 0;

    /**
     * Attach source.
     * Given a vector of abstract source objects, this function will query each for its spatial location. After
     * performing a convex hull test, it will perform a quick inverse problem to determine the position of any sources
     * within each element in reference coordinates. These reference coordinates are then saved in the source object.
     * References to any sources which lie within the element are saved in the mSources vector.
     * @param [in] sources A vector of all the sources defined for a simulation run.
     */
    virtual void attachSource(std::vector<Source*> sources) = 0;

    /**
     * Simple function to set the (remembered) element number.
     */
    virtual void SetLocalElementNumber(int element_number) { mElementNumber = element_number; }
    
    /**
     * Sums a field into the mesh (global dofs) owned by the current processor.
     * This function sets up and calls PLEX's DMVecSetClosure for a given element. Remapping is handled implicitly.
     * @param mesh [in] A reference to the mesh to which this element belongs.
     * @param field [in] The values of the field on the element, in Salvus ordering.
     * @param name [in] The name of the global fields where the field will be summed.
     * TODO: Make this function check if the field is valid?
     */
    virtual void checkInFieldElement(Mesh *mesh, const Eigen::VectorXd &field, const std::string name) = 0;
    
    /**
     * Sets a field into the mesh (global dofs) owned by the current processor.
     * This function sets up and calls PLEX's DMVecSetClosure for a given element. Remapping is handled implicitly.
     * @param mesh [in] A reference to the mesh to which this element belongs.
     * @param field [in] The values of the field on the element, in Salvus ordering.
     * @param name [in] The name of the global fields where the field will be summed.
     * TODO: Make this function check if the field is valid?
     */
    virtual void setFieldElement(Mesh *mesh, const Eigen::VectorXd &field, const std::string name) = 0;
    
    /**
     * Queries the mesh for, and returns, a field.
     * This function returns a Matrix (or Vector -- a Eigen::Vector is just a special case of a matrix) from the
     * global dofs owned by the current processor. The returned field will be in Salvus ordering -- remapping is
     * done implicitly. If the field is multi dimensional, the dimensions will be ordered as rows (i.e.
     * row(0)->x, row(1)->z).
     * @param [in] mesh Pointer to the mesh representing the current element.
     * @param [in] name Name of field to check out.
     */
    virtual Eigen::VectorXd checkOutFieldElement(Mesh *mesh, const std::string name) = 0;

    /**
     * Builds nodal coordinates (x,z) on all mesh degrees of freedom.
     * @param mesh [in] The mesh.
     */
    virtual std::tuple<Eigen::VectorXd,Eigen::VectorXd> buildNodalPoints(Mesh* mesh) = 0;

    /**
     * Figure out which dofs (if any) are on the boundary.
     */
    virtual void setBoundaryConditions(Mesh *mesh) = 0;


    // Attribute gets.
    virtual int Number() const { return mElementNumber; }
    virtual int NumberDofEdge() const { return mNumberDofEdge; }
    virtual int NumberDofFace() const { return mNumberDofFace; }
    virtual int NumberDofVertex() const { return mNumberDofVertex; }        
    virtual int NumberIntegrationPoints() const { return mNumberIntegrationPoints; }

    static Eigen::VectorXi ElementClosure() { return mClosureMapping; }
    
    int NumberDimensions() const { return mNumberDimensions; }
    
    virtual std::map<std::string,int> Boundaries() const { return mBoundaries; }

    virtual Eigen::MatrixXd GetVertexCoordinates() { return mVertexCoordinates; }


    // Abstract methods to be implemented by elements at the physics level (acoustic/elastic)
    virtual Eigen::MatrixXd computeSourceTerm(double time) = 0;
    virtual void assembleElementMassMatrix(Mesh *mesh) = 0;
    virtual void interpolateMaterialProperties(ExodusModel *model) = 0;
    virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement) = 0;

    virtual std::vector<std::string> PullElementalFields() const = 0;
    virtual std::vector<std::string> PushElementalFields() const = 0;

};
