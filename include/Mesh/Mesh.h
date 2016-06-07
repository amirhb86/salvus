#pragma once

// stl.
#include <set>
#include <map>
#include <iosfwd>
#include <string>
#include <vector>
#include <assert.h>
#include <iostream>

// 3rd party.
#include <petsc.h>
#include <Eigen/Dense>

// forward decl.
class Options;
class ExodusModel;

/**
 * Struct holding the vectors representing the global DOFs.
 * The global PETSc vector is defined across processors,
 * and knows which of its entries are on processor boundaries. The local PETSc vector exists only on the
 * current processor.
 */
struct vec_struct {
  std::string name;
  /** < Field name (i.e. displacement_x) */
  Vec glb;
  /** < Global PETSc vector */
  Vec loc;            /** < Local PETSc vector */
};

class Mesh {

  /** Keeps track of all the fields defined in the mesh. **/
  std::set<std::string> mMeshFields;

  /** Keeps track of the primary field on each element. **/
  std::vector<std::set<std::string>> mElmFields;

  /** Keeps track of any (possible) coupling fields on each element. **/
  std::map<PetscInt,std::set<std::string>> mCouplingFields;

  /** Keeps (sparse) track of any specific couplings. **/
  std::map<PetscInt,std::vector<std::tuple<PetscInt,std::vector<std::string>>>> mCpl;

  int mNumberElementsLocal;
  /** < Num of elements on this processor. */
  int mNumberDimensions;
  /** < Num of dimensions of the mesh. */
  int mNumberSideSets;
  /** < Num of flagged boundaries. */

  std::string mExodusFileName;
  /** < Exodus file from which this mesh (skeleton) was defined. */

  DM mDistributedMesh;
  /** < PETSc distributed mesh defining parallel element layout. */
  PetscSection mMeshSection;      /** < Mesh section describing location of the integration points. In the future we
                                        * may have many of these per mesh. */
  int int_tstep;


 protected:

  std::vector<std::string> mGlobalFields;
  /** < List of field names on global dof ("u","v",etc) */
  std::map<std::string, vec_struct> mFields;
  /** < Dictionary holding the fields on the global dof. */
  PetscViewer mViewer;
  /** < Holds information used to dump field values to disk. */

  /** The CFL constant for the time-stepping scheme (Newmark 2nd-order sets 1.0)*/
  double mCFL = 1.0;

  std::map<PetscInt, std::string> mBoundaryIds;
  /** < mapping between boundary id
                                                       and its associated name (e.g.,
                                                       "absorbing boundary" in the
                                                       petsc side set collection) */


  std::map<std::string, std::map<int, std::vector<int>>>
      mBoundaryElementFaces;  /** < list of elements on a boundary and the corresponding
                                    boundary faces. Each boundary has its own list of
                                    elements, and each element has its own list of
                                    faces. */


 public:

  /**
   * Factory which returns a `Mesh` based on user defined options.
   * The `Mesh` is really a collection of the global dofs, and the fields defined on these dofs will depend on both
   * the physics of the system under consideration, and the method of time-stepping chosen.
   * @return Some derived mesh class.
   */
  static Mesh *factory(Options options);

  /**
   * Given an existing vector of continuous fields, append a new set of fields based on a
   * type of physics.
   * @param [in] fields A pre-existing vector containing field names.
   * @param [in] physics A string defining the physics you would like to add.
   * @return An extended vector containing new field names.
   */
  static std::vector<std::string> appendPhysicalFields(const std::vector<std::string>& fields,
                                                       const std::string& physics);

  virtual ~Mesh() { DMDestroy(&mDistributedMesh); }

  /**
   * Reads an exodus mesh from a file defined in options.
   * By the time this method is finished, the mesh has been
   * read, and parallelized across processors (via a call to Chaco).
   * @param [in] options The master options struct. TODO: Move the gobbling of options to the constructor.
   */
  void read(Options options);

  /**
   * Alternate version of read, which creates a mesh given a matrix of verts and cells.
   */
  void read(int dim, int numCells, int numVerts, int numVertsPerElem,
            Eigen::MatrixXi cells, Eigen::MatrixXd coords);
  
  /**
   * Sets up the dofs across elements and processor boundaries.
   * Specifically, this function defines a `DMPlex section` spread across all elements. It is this section that
   * defines the integration points, and carries information about which of these points are shared across processor
   * boundaries. This function takes a model object because it needs to know about the physics attached
   * to each element. It uses the physics to set an appropriate amount of components belonging to each DOF.
   * @param [in] number_dof_vertex Number of dofs per 0-d mesh component (vertex). 1 for the standard GLL basis.
   * @param [in] number_dof_edge Number of dofs per 1-d mesh component (edge). order-1 for the standard GLL basis.
   * @param [in] number_dof_face Number of dofs per 2-d mesh component (face). (order-1)^2 for the standard GLL basis.
   * @param [in] number_dof_volume Num of dofs per 3-d mesh component (volume). Something something for the
   * standard GLL basis.
   */
  PetscErrorCode setupGlobalDof(int number_dof_vertex,
                                int number_dof_edge,
                                int number_dof_face,
                                int number_dof_volume,
                                int number_dimensions,
                                ExodusModel *model);

  /**
   * Registers both the global (across parallel partition) and local (on a single parallel partition) vectors for a
   * given name. Vector information is stored in an `std::map` under mFields, with `name` as the key.
   * @param name Name of field to register.
   */
  void registerFieldVectors(const std::string &name);

  /**
   * Begins (and ends) the gloabl -> local MPI sends for a given field name. By the time this function returns, you
   * can be confident that the MPI sends have been completed (i.e. it is blocking), and will be available to each
   * element.
   * @param name Name of field to checkout.
   */
  void checkOutField(const std::string &name);

  /**
   * Does the local -> global MPI sends for a given field name. The
   * send is performed with an sum, i.e.  the value of field on
   * coincident GLL points are properly summed together.
   * @param name Name of field to assemble on global dofs.
   */
  void assembleLocalFieldToGlobal(const std::string &name);

  /**
   * Begins the local -> global MPI sends for a given field name. The
   * send is performed with an sum, i.e.  the value of field on
   * coincident GLL points are properly summed together. Note that
   * this function is non-blocking!! This MUST be paired with an
   * equivalent call to `assembleLocalFieldToGlobalEnd`.
   * @param name Name of field to assemble on global dofs.
   */
  void assembleLocalFieldToGlobalBegin(const std::string &name);

  /**
   * Makes a processor local array "global". As a "set", does not incurr any communication.
   * @param name Name of field to assemble on global dofs.
   */
  void setLocalFieldToGlobal(const std::string &name);

  /**
   * Finishes the local -> global MPI sends for a given field name. The
   * send is performed with an sum, i.e.  the value of field on
   * coincident GLL points are properly summed together. This MUST be paired with an
   * equivalent call to `assembleLocalFieldToGlobalBegin`.
   * @param name Name of field to assemble on global dofs.
   */
  void assembleLocalFieldToGlobalEnd(const std::string &name);

  /**
   * Begins the local -> global MPI sends for a given field name. The send is performed with an implied sum, i.e.
   * the value of field on coincident GLL points are properly summed together. Note that this function is
   * non-blocking!! This MUST be paired with an equivalent call to checkInFieldEnd.
   * @param name Name of field to checkin.
   */
  void checkInFieldBegin(const std::string &name);

  /**
   * Ends the local -> global MPI send for a given field name. This function should come after an equivalent
   * checkInFieldBegin. This method is blocking, so you can be confident that when it returns the desired field has
   * been scattered and summed into the global (parallel) degrees of freedom.
   * @param name Name of field to checkin.
   */
  void checkInFieldEnd(const std::string &name);

  /**
   * Returns an ordered vector of a field (i.e. x-displacement) on a
   * PETSc point (e.g., element,face,vertex,etc), via a call to
   * `DMPlexVecGetClosure`.
   * @param [in] point The PETSc point
   * @param [in] name Name of field.
   * @ return The field on the point
   */
  Eigen::VectorXd getFieldOnPoint(int point,std::string name);

  
  /**
   * Returns an ordered vector of a field (i.e. x-displacement) on an element, via a call to DMPlexVecGetClosure.
   * Note that a vector containing the closure mapping must also be passed in -- this should change in the future.
   * @param [in] name Name of field.
   * @param [in] element_number Element number (on the local processor) for which the field is requested.
   * @param [in] closure A vector of size mNumIntPnt, which specifies the mapping from the Plex element
   * closure to the desired gll point ordering.
   * @ return The ordered field on an element.
   */
  Eigen::VectorXd getFieldOnElement(const std::string &name, const int &element_number,
                                    const Eigen::VectorXi &closure);

  /**
   * Returns an ordered vector of a field (i.e. x-displacement) on a face, via a call to
   * DMPlexVecGetClosure.  Note that a vector containing the closure mapping (for the face) must
   * also be passed in -- this should change in the future.
   * @param [in] name Name of field.
   * @param [in] element_number Element number (on the local processor) for which the field is requested.
   * @ return The ordered field on an element.
   */
  Eigen::VectorXd getFieldOnFace(const std::string &name, const int &face_number);

  /**
   * Sets a field from a face into the degrees of freedom owned by the local processor, via a call to
   * DMPlexVecSetClosure. Shared DOF should be identical, and are thus set from an arbitrary element.
   * @param [in] field_name Name of field.
   * @param [in] face_number Face number (on the local processor) for which the field is requested.
   * @param [in] closure A vector of size mNumIntPnt, which specifies the mapping from the Plex face
   * closure to the desired gll point ordering.
   * @param [in] field The element-ordered field (i.e. x-displacement) to insert into the mesh.
   */
  void setFieldFromFace(const std::string &name, const int face_number, const Eigen::VectorXd &field);

  /**
   * Adds a field from a face into the degrees of freedom owned by the local processor, via a call to
   * DMPlexVecSetClosure. Shared DOF are "assembled" (added), and are thus a sum from each arbitrary element.
   * @param [in] field_name Name of field.
   * @param [in] face_number Face number (on the local processor) for which the field is requested.
   * @param [in] closure A vector of size mNumIntPnt, which specifies the mapping from the Plex face
   * closure to the desired gll point ordering.
   * @param [in] field The element-ordered field (i.e. x-displacement) to insert into the mesh.
   */
  void addFieldFromFace(const std::string &name, const int face_number, const Eigen::VectorXd &field);

  /**
   * Sets a field from an element into the degrees of freedom owned by the local processor, via a call to
   * DMPlexVecSetClosure. Shared DOF should be identical, and are thus set from an arbitrary element.
   * @param [in] field_name Name of field.
   * @param [in] element_number Element number (on the local processor) for which the field is requested.
   * @param [in] closure A vector of size mNumIntPnt, which specifies the mapping from the Plex element
   * closure to the desired gll point ordering.
   * @param [in] field The element-ordered field (i.e. x-displacement) to sum into the mesh.
   */
  void setFieldFromElement(const std::string &name, const int element_number,
                           const Eigen::VectorXi &closure, const Eigen::VectorXd &field);

  /**
   * Sums a field from an element into the degrees of freedom owned by the local processor, via a call to
   * DMPlexVecSetClosure. Shared DOF are summed (i.e., assembled).
   * @param [in] field_name Name of field.
   * @param [in] element_number Element number (on the local processor) for which the field is requested.
   * @param [in] closure A vector of size mNumIntPnt, which specifies the mapping from the Plex element
   * closure to the desired gll point ordering.
   * @param [in] field The element-ordered field (i.e. x-displacement) to sum into the mesh.
   */
  void addFieldFromElement(const std::string &name, const int element_number,
                           const Eigen::VectorXi &closure, const Eigen::VectorXd &field);

  /**
   * Number of elements owned by the current processors.
   * @ return Num of Elements.
   */
  inline PetscInt NumberElementsLocal() { return mNumberElementsLocal; }

  /**
   * Defines the field to visualize, and the name of the file to dump to.
   * @param [in] movie_filename Where to save the file.
   * TODO: Should add field name here, and interval? Do we need to support other dump options than HDF5?
   */
  void setUpMovie(const std::string &movie_filename);

  /**
   * Commits a movie frame to disk, using the interface defined in setUpMovie.
   */
  void saveFrame(std::string name, PetscInt timestep);

  /**
   * Needs to be called at the end of a time loop if a movie is desired.
   * Cleans up the open HDF5 file, and safely shuts down the output stream.
   */
  void finalizeMovie();

  /**
   * Set a desired vec_struct to zero. Usefull for those terms which we sum into (i.e. element forcing).
   */
  void zeroFields(const std::string &name);

  /**
   * Read and setup boundaries (exodus side sets). Gets list of faces from Petsc and builds
   * corresponding list of elements, which lay on each boundary.
   */
  int setupBoundaries(Options options);

  /**
   * Reads boundaries (exodus side sets) from mesh file and builds appropriate boundary name to
   * petsc id mapping.
   */
  int readBoundaryNames(Options options);

  /**
   * Virtual function which implements the time stepping.
   * This will change based on physics, time stepping routine, and dimension.
   */
  virtual void advanceField(double dt) = 0;

  /**
   * Virtual function to apply the inverse of a mass matrix.
   * This is usually called at the beginning of a time loop. Again, it will depend on the specific physics and
   * time scheme used. When this function is complete, the global fields should be ready to take a time-step (at
   * least that's the case for explicit schemes).
   */
  virtual void applyInverseMassMatrix() = 0;

  /**
   * Return list of all the fields required on the global degrees of freedom.
   */
  std::vector<std::string> GlobalFields() const { return mGlobalFields; }

  /**
   * Add field to list of global degrees of freedom. Has to be done
   * before mesh initializes fields to global DOF. (via registerFieldVectors)
   */
  void AddToGlobalFields(std::string fieldname) { mGlobalFields.push_back(fieldname); }

  PetscInt GetNeighbouringElement(const PetscInt interface, const PetscInt this_elm) const;

  /**
   * Get the transitive closure of a coordinate section for a mesh point.
   */
  Eigen::MatrixXd getElementCoordinateClosure(PetscInt elem_num);

  int numFieldPerPhysics(std::string physics);

  inline std::vector<std::string> ElementFields(const PetscInt num) {
    return std::vector<std::string> (mElmFields[num].begin(), mElmFields[num].end());
  }

  inline DM &DistributedMesh() { return mDistributedMesh; }
  inline PetscSection &MeshSection() { return mMeshSection; }
  virtual std::map<PetscInt, std::string> &BoundaryIds() { return mBoundaryIds; }

  inline int NumberSideSets() { return mNumberSideSets; }
  inline int NumberDimensions() { return mNumberDimensions; }

  inline std::map<std::string, std::map<int, std::vector<int>>>
  BoundaryElementFaces() { return mBoundaryElementFaces; }
  std::set<std::string> AllFields() const { return mMeshFields; }

  std::vector<std::tuple<PetscInt,std::vector<std::string>>> CouplingFields(const PetscInt elm);
  std::vector<std::string> TotalCouplingFields(const PetscInt elm);
  std::vector<PetscInt> EdgeNumbers(const PetscInt elm);

  inline double CFL() { return mCFL; }
  
};
