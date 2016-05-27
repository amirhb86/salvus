#pragma once

// stl.
#include <vector>
#include <assert.h>

// 3rd party.
#include <mpi.h>
#include <petsc.h>
#include <Eigen/Dense>

extern "C" {
#include "Utilities/kdtree.h"
#include "exodusII.h"
};

// forward decl.
class Options;
class Utilities;

class ExodusModel {

  int mExodusId;
  int mNumberDimension;
  int mNumberVertices;
  int mNumberElements;
  int mNumberElementBlocks;
  int mNumberNodeSets;
  int mNumberSideSets;
  int mNumberNodalVariables;
  int mNumberElementVertex;

  std::vector<int> mElementBlockIds;
  std::vector<int> mVerticesPerElementPerBlock;
  
    // New variables for pymesher.
    int mNumberElementalVariables;
    std::vector<std::string> mElementalVariableNames;
    std::vector<double> mElementalVariables;
    kdtree *mElementalKdTree;
    std::vector<int> mElementalKdTreeData;

  char mExodusTitle[MAX_LINE_LENGTH + 1];

  float mExodusVersion;

  kdtree *mNodalKdTree;

  std::string mExodusFileName;

  std::vector<int> mElementConnectivity;
  std::vector<int> mNodalKdTreeData;
  std::vector<double> mNodalX;
  std::vector<double> mNodalY;
  std::vector<double> mNodalZ;
  std::vector<double> mNodalVariables;

  std::vector<std::string> mNodalVariableNames;
  std::vector<std::string> mSideSetNames;

  void getInitialization();
  void readCoordinates();
  void readNodalVariables();
  void readElementalVariables();
  void createNodalKdTree();
  void createElementalKdTree();
  void readConnectivity();
  void exodusError(const int retval, std::string func_name);

  /**
   * Query the sidesets from the exodus file. Save the names present into the mSideSetNames vector.
   */
  void readSideSets();


 public:

  ExodusModel(Options options);
  ~ExodusModel() {
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (!rank) { ex_close(mExodusId); }
  }
  void initializeParallel();
  PetscReal getMaterialParameterAtPoint(const std::vector<double> point,
                                        const std::string parameter_name);

  /**
   * Returns a parameter at a specific vertex, following the ordering in the reference element. For example,
   * in a 2D quad, we have P_0, P_1, P_2, P_3. See the reference element ordering for further details.
   * @param [in] elem_center Vector containing an ordered pair of doubles representing the element center.
   * @param [in] parameter_name The base parameter name.
   * @param [in] vertex_num The vertex for which the parameter is desired.
   * @return The value of the material parameter at the specified vertex.
   */
  double getElementalMaterialParameterAtVertex(const Eigen::VectorXd &elem_center,
                                               const std::string &parameter_name,
                                               const int vertex_num) const;

  /**
   * Returns a string specifying which type of material the element is.
   * I.e. will return "ACOUSTIC" for acoustic, "ELASTIC" for elastic.
   */
  std::string getElementType(const Eigen::VectorXd &elem_center);

  int NumberElementVertices() const { return mNumberElementVertex; }

};

