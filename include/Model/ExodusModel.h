#pragma once

// stl.
#include <vector>
#include <assert.h>
#include <memory>

// 3rd party.
#include <mpi.h>
#include <petsc.h>
#include <Eigen/Dense>
#include <Utilities/Types.h>

extern "C" {
#include "Utilities/kdtree.h"
#include "exodusII.h"
};

// forward decl.
class Options;
class Utilities;

class ExodusModel {

 protected:

  /* Title stored in file. */
  char mExodusTitle[MAX_LINE_LENGTH + 1];
  std::string mExodusFileName;

  /* Integer descriptions. */
  PetscInt mExodusId;
  PetscInt mNumberNodeSets;
  PetscInt mNumberSideSets;
  PetscInt mNumberVertices;
  PetscInt mNumberElements;
  PetscInt mNumberDimension;
  PetscInt mNumberElementBlocks;
  PetscInt mNumberNodalVariables;
  PetscInt mNumberElementalVariables;
  PetscInt mNumberGlobalVariables;
  PetscInt mNumberInfo;

  /* Version number. */
  float mExodusVersion;

  /* Functionality for multiple blocks. */
  std::vector<PetscInt> mElementBlockIds;
  std::vector<PetscInt> mElementConnectivity;
  std::vector<PetscInt> mVerticesPerElementPerBlock;

  /* kDtree based on either element centres, or vertices. */
  kdtree *mNodalKdTree;
  kdtree *mElementalKdTree;
  std::vector<PetscInt> mNodalKdTreeData;
  std::vector<PetscInt> mElementalKdTreeData;

  /* Vector to hold variables on each element. */
  std::vector<std::string> mElementalVariableNames;
  std::vector<PetscReal> mElementalVariables;

  /* Vector to hold variables defined only at nodes. */
  std::vector<PetscReal> mNodalVariables;
  std::vector<std::string> mNodalVariableNames;

  /* Vector to hold global mesh variables. */
  std::vector<PetscReal> mGlobalVariables;
  std::vector<std::string> mGlobalVariableNames;

  /* Info */
  std::vector<std::string> mInfo;

  /* Nodal locations. */
  std::vector<PetscReal> mNodalX;
  std::vector<PetscReal> mNodalY;
  std::vector<PetscReal> mNodalZ;

  /* Side sets define edge boundary conditions. */
  std::vector<std::string> mSideSetNames;

  /** Read (dimension specific) coordinate values. */
  void readCoordinates();
  /** Read mesh connectivity. */
  void readConnectivity();
  /** Get parameters from Exodus file. */
  void getInitialization();
  /** Create kDTree based on element vertices. */
  void createNodalKdTree();
  /** Read variables defined at nodal locations. */
  void readNodalVariables();
  /** Read global variables. */
  void readGlobalVariables();
  /** Create kDtree based on element centers. */
  void createElementalKdTree();
  /** Read variables defined elementwise (i.e. VP_0, VP_1, ... ,VP_n). */
  void readElementalVariables();
  /** Query the sidesets from the exodus file. Save the names present into the mSideSetNames vector. */
  void readSideSets();
  /** Read Exodus info field. */
  void readInfo();

  /**
   * Throw an error and specify which function has errored.
   * @param [in] retval Value thrown by exodus library function.
   * @param [in] func_name Name of function called.
   */
  void exodusError(const int retval, std::string func_name);

 public:

  /** Constructor (without specifying the filename yet). */
  ExodusModel();
  /** Constructor (sets filename). */
  ExodusModel(std::unique_ptr<Options> const &options);

  /** Destructor (closes exodus file). */
  ~ExodusModel();

  /** Set filename to read the model from. */
  void setExodusFilename(const std::string filename);

  /**
   * Reads in mesh on rank 0, and populates all necessary model quantities (parameters, etc.).
   * Broadcasts the necessary parameters to all processors.
   */
  void read();

  /**
   * Writes out mesh on rank 0, including all necessary model quantities (parameters, etc.).
   */
  void write(const std::string filename);

  /**
   * Returns the closest parameter to a point (i.e. element node).
   * @param [in] point Point which to search for.
   * @param [in] parameter_name Name of material parameter.
   * @return The value of the closest material parameter defined at a point.
   */
  PetscReal getNodalParameterAtNode(const Eigen::Ref<const RealVec>& point,
                                    const std::string parameter_name);

  /**
   * Returns a parameter at a specific vertex, following the ordering in the reference element. For example,
   * in a 2D quad, we have P_0, P_1, P_2, P_3. See the reference element ordering for further details.
   * @param [in] elem_center Vector containing an ordered pair of doubles representing the element center.
   * @param [in] parameter_name The base parameter name.
   * @param [in] vertex_num The vertex for which the parameter is desired.
   * @return The value of the material parameter at the specified vertex.
   */
  PetscScalar getElementalMaterialParameterAtVertex(const Eigen::Ref<const RealVec>& elem_center,
                                                    std::string parameter_name,
                                                    const PetscInt vertex_num) const;

  /**
   * Returns a string specifying which type of material the element is.
   * I.e. will return "ACOUSTIC" for acoustic, "ELASTIC" for elastic.
   */
  std::string getElementType(const Eigen::VectorXd &elem_center);

  /**
   * Returns the a string specifying the side set name (x0, x1, ...).
   * @param side_set_num Exodus side set identifier [0, nSideSet].
   * @return Side set name.
   */
  std::string SideSetName(const PetscInt side_set_num);


};

