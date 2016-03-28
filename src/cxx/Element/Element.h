#pragma once

#include <Eigen/Dense>
#include <Model/ExodusModel.h>
#include <Source/Source.h>
#include <Mesh/Mesh.h>
#include <mpi.h>

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
    


    /**********************************************************************************
     * OBJECT MEMBERS. THESE VARIABLES AND FUNCTIONS SHOULD APPLY TO SPECIFIC ELEMENTS.
     ***********************************************************************************/

    int mElementNumber; /** Element number on the local processor. */

    std::vector<Source*> mSources;  /** Vector of abstract sources belonging (spatial) to the element */

    Eigen::VectorXd mMassMatrix;    /** Elemental mass matrix */
    
    Eigen::MatrixXd mVertexCoordinates;   /** Vertex coordinates ordered as above. row(0)->x, row(1)->z */

    Eigen::Matrix<double,2,1> mElementCenter; /** (x, z) location of element center */

    bool mOnBoundary;    /** < Whether or not the current element has a special boundary condition */

    std::map<std::string,std::vector<int>> mBoundaries;  /** < Map relating the type of boundary to the Petsc edge number */

public:

    /**
     * Factory return the proper element physics based on the command line options.
     * @return Some derived element class.
     */
    static Element2D *factory(Options options);
    
    /**
     * Copy constructor.
     * Returns a copy. Use this once the reference element is set up via the constructor, to allocate space for
     * all the unique elements on a processor.
     */
    virtual Element2D * clone() const = 0;
    
    /**
     * Utility function to integrate a field over the element. This could probably be made static, but for now I'm
     * just using it to check some values.
     * @param [in] field The field which to integrate, defined on each of the gll points.
     * @returns The scalar value of the field integrated over the element.
     */
    virtual double integrateField(const Eigen::VectorXd &field) = 0;

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
    virtual void checkInFieldElement(Mesh *mesh, const Eigen::VectorXd &field, const std::string name);
    
    /**
     * Sets a field into the mesh (global dofs) owned by the current processor.
     * This function sets up and calls PLEX's DMVecSetClosure for a given element. Remapping is handled implicitly.
     * @param mesh [in] A reference to the mesh to which this element belongs.
     * @param field [in] The values of the field on the element, in Salvus ordering.
     * @param name [in] The name of the global fields where the field will be summed.
     * TODO: Make this function check if the field is valid?
     */
    virtual void setFieldElement(Mesh *mesh, const Eigen::VectorXd &field, const std::string name);
    
    /**
     * Queries the mesh for, and returns, a field.
     * This function returns a Matrix (or Vector -- a Eigen::Vector is just a special case of a matrix) from the
     * global dofs owned by the current processor. The returned field will be in Salvus ordering -- remapping is
     * done implicitly. If the field is multi dimensional, the dimensions will be ordered as rows (i.e.
     * row(0)->x, row(1)->z).
     * @param [in] mesh Pointer to the mesh representing the current element.
     * @param [in] name Name of field to check out.
     */
    virtual Eigen::VectorXd checkOutFieldElement(Mesh *mesh, const std::string name);

    /**
     * Builds nodal coordinates (x,(y),z) on all mesh degrees of freedom.
     * @param mesh [in] The mesh.
     */
    virtual std::tuple<Eigen::VectorXd,Eigen::VectorXd> buildNodalPoints(Mesh* mesh) = 0;

    // Attribute gets.
    virtual int Number() const { return mElementNumber; }
    virtual int NumberDofEdge() const { return mNumberDofEdge; }
    virtual int NumberDofFace() const { return mNumberDofFace; }
    virtual int NumberDofVertex() const { return mNumberDofVertex; }        
    virtual int NumberIntegrationPoints() const { return mNumberIntegrationPoints; }
    
    
    virtual Eigen::VectorXi ElementClosure() { return mClosureMapping; }
    
    int NumberDimensions() const { return mNumberDimensions; }
    
    virtual std::map<std::string,std::vector<int>> Boundaries() const { return mBoundaries; }

    virtual Eigen::MatrixXd GetVertexCoordinates() { return mVertexCoordinates; }


    // Abstract methods to be implemented by elements at the physics level (acoustic/elastic)
    virtual Eigen::MatrixXd computeSourceTerm(double time) = 0;
    virtual void assembleElementMassMatrix(Mesh *mesh) = 0;
    virtual void interpolateMaterialProperties(ExodusModel *model) = 0;
    virtual Eigen::MatrixXd computeStiffnessTerm(const Eigen::MatrixXd &displacement) = 0;
    virtual void prepareStiffness() = 0;

    inline bool OnBoundary() { return mOnBoundary; }
    
    /**
     * Figure out which dofs (if any) are on the boundary.
     */
    virtual void setBoundaryConditions(Mesh *mesh);

    /**
     * Apply boundary condition (dirichlet, absorbing, fluid-elastic coupling, etc)
     */ 
    virtual void applyBoundaryConditions(Mesh *mesh,
                                         Options &options,
                                         std::string fieldname);
    
    
    virtual std::vector<std::string> PullElementalFields() const = 0;
    virtual std::vector<std::string> PushElementalFields() const = 0;

    // Testing, which is usually implemented via an initial condition
    virtual void setupTest(Mesh* mesh, Options options)  {
//        printf("ERROR: No test implemented\n");
//        MPI::COMM_WORLD.Abort(-1);
    }
    virtual double checkTest(Mesh* mesh, Options options, const Eigen::MatrixXd &displacement, double time) {
//        printf("ERROR: No test implemented\n");
//        MPI::COMM_WORLD.Abort(-1);
//        return -1;
    }
    
    
};
