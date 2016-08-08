#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Utilities/Utilities.h>
#include <Utilities/Logging.h>
#include <Utilities/Types.h>
#include <algorithm>

ExodusModel::ExodusModel(std::unique_ptr<Options> const &options) {
  mNodalKdTree = NULL;
  mElementalKdTree = NULL;
  mExodusId = 0;
  mExodusFileName = options->ModelFile();
}

ExodusModel::~ExodusModel() {
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank && mExodusId) { ex_close(mExodusId); }
  if (mElementalKdTree) { kd_free(mElementalKdTree); }
  if (mNodalKdTree) { kd_free(mNodalKdTree); }
}

void ExodusModel::read() {

  /* Read the model from rank 0. */
  int root = 0;
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank)
  {
    getInitialization();
    readConnectivity();
    readCoordinates();
    readElementalVariables();
    readNodalVariables();
    readSideSets();
  }

  /* Broadcast all scalars. */
  mNumberElements = utilities::broadcastNumberFromRank(mNumberElements, root);
  mNumberNodeSets = utilities::broadcastNumberFromRank(mNumberNodeSets, root);
  mNumberSideSets = utilities::broadcastNumberFromRank(mNumberSideSets, root);
  mNumberVertices = utilities::broadcastNumberFromRank(mNumberVertices, root);
  mNumberDimension = utilities::broadcastNumberFromRank(mNumberDimension, root);
  mNumberElementBlocks = utilities::broadcastNumberFromRank(mNumberElementBlocks, root);

  /* Broadcast all vectors. */
  mNodalX = utilities::broadcastNumberVecFromRank(mNodalX, root);
  mNodalY = utilities::broadcastNumberVecFromRank(mNodalY, root);
  mSideSetNames = utilities::broadcastStringVecFromRank(mSideSetNames, root);
  mElementalVariables = utilities::broadcastNumberVecFromRank(mElementalVariables, root);
  mElementConnectivity = utilities::broadcastNumberVecFromRank(mElementConnectivity, root);
  mElementalVariableNames = utilities::broadcastStringVecFromRank(mElementalVariableNames, root);
  mVerticesPerElementPerBlock = utilities::broadcastNumberVecFromRank(mVerticesPerElementPerBlock, root);

  /* Broadcast dimension specific components. */
  if (mNumberDimension > 2) {
    mNodalZ = utilities::broadcastNumberVecFromRank(mNodalZ, root);
  }

  /* Element tree labels elements by centroid. Nodal tree labels elements by vertex. */
  createElementalKdTree();
  createNodalKdTree();

}

void ExodusModel::exodusError(const int retval, std::string func_name) {

  if (retval < 0) {
    throw std::runtime_error(
        "Error in exodus function: " + func_name + " with retval " + std::to_string((long long) retval));
  } else if (retval > 0) {
    throw salvus_warning(
        "Warning from exodus function: " + func_name);
  }

}

void ExodusModel::readCoordinates() {

  mNodalX.resize(mNumberVertices);
  mNodalY.resize(mNumberVertices);
  if (mNumberDimension == 3) { mNodalZ.resize(mNumberVertices); }
  if (mNumberDimension == 2) {
    exodusError(ex_get_coord(
        mExodusId, mNodalX.data(), mNodalY.data(), NULL),
                "ex_get_coord");
  } else {
    exodusError(ex_get_coord(
        mExodusId, mNodalX.data(), mNodalY.data(), mNodalZ.data()),
                "ex_get_coord");
  }

}

void ExodusModel::getInitialization() {

  int io_ws = 0;
  int comp_ws = 8;
  mExodusId = ex_open(mExodusFileName.c_str(), EX_READ, &comp_ws, &io_ws, &mExodusVersion);
  if (mExodusId < 0) { throw std::runtime_error("Error opening exodus model file."); }

  exodusError(ex_get_init(
      mExodusId, mExodusTitle, &mNumberDimension, &mNumberVertices, &mNumberElements,
      &mNumberElementBlocks, &mNumberNodeSets, &mNumberSideSets),
              "ex_get_init");
  mElementBlockIds.resize(mNumberElementBlocks);
  mVerticesPerElementPerBlock.resize(mNumberElementBlocks);
  exodusError(ex_get_elem_blk_ids(mExodusId, mElementBlockIds.data()),
              "ex_get_elem_blk_ids");
  for (int i = 0; i < mNumberElementBlocks; i++) {
    char elem_type[256];
    int num_elem_this_blk;
    int num_nodes_per_elem;
    int num_attr;
    exodusError(ex_get_elem_block(mExodusId,
                                  mElementBlockIds[i],
                                  elem_type, &num_elem_this_blk,
                                  &num_nodes_per_elem,
                                  &num_attr
    ), "ex_get_elem_block");
    mVerticesPerElementPerBlock[i] = num_nodes_per_elem;
  }

}

void ExodusModel::createNodalKdTree() {

  // Need to keep the index arrays allocated.
  mNodalKdTreeData.resize(mNodalX.size());

  // Create.
  if (mNumberDimension == 2) {
    mNodalKdTree = kd_create(2);
    for (auto i = 0; i < mNumberVertices; i++) {
      mNodalKdTreeData[i] = i;
      kd_insert(mNodalKdTree, std::vector<double> {mNodalX[i], mNodalY[i]}.data(), &mNodalKdTreeData[i]);
    }
  } else {
    mNodalKdTree = kd_create(3);
    for (auto i = 0; i < mNumberVertices; i++) {
      mNodalKdTreeData[i] = i;
      kd_insert(mNodalKdTree, std::vector<double> {mNodalX[i], mNodalY[i], mNodalZ[i]}.data(),
                &mNodalKdTreeData[i]);
    }
  }

}


void ExodusModel::createElementalKdTree() {

  // Need to keep the index arrays allocated.
  mElementalKdTreeData.resize(mNumberElements);

  // TODO: Multiple blocks?
  int num_vertex_per_elem = mVerticesPerElementPerBlock[0]; // first block

  if (mNumberDimension == 2) {
    mElementalKdTree = kd_create(2);
    for (auto i = 0; i < mNumberElements; i++) {

      // Calculate element centers.
      double avg_x = 0, avg_y = 0;
      for (auto j = num_vertex_per_elem * i; j < num_vertex_per_elem * (i + 1); j++) {
        avg_x += mNodalX[mElementConnectivity[j] - 1];
        avg_y += mNodalY[mElementConnectivity[j] - 1];
      }

      // Insert centers into tree.
      mElementalKdTreeData[i] = i;
      kd_insert(mElementalKdTree, std::vector<double> {
                    avg_x / (double) num_vertex_per_elem, avg_y / (double) num_vertex_per_elem}.data(),
                &mElementalKdTreeData[i]);
    }

  } else if (mNumberDimension == 3) {
    mElementalKdTree = kd_create(3);
    for (auto i = 0; i < mNumberElements; i++) {

      // Calculate element centers.
      double avg_x = 0, avg_y = 0, avg_z = 0;
      for (auto j = num_vertex_per_elem * i; j < num_vertex_per_elem * (i + 1); j++) {
        avg_x += mNodalX[mElementConnectivity[j] - 1];
        avg_y += mNodalY[mElementConnectivity[j] - 1];
        avg_z += mNodalZ[mElementConnectivity[j] - 1];
      }

      // Insert centers into tree.
      mElementalKdTreeData[i] = i;
      kd_insert(mElementalKdTree,
                std::vector<double> {
                    avg_x / (double) num_vertex_per_elem,
                    avg_y / (double) num_vertex_per_elem,
                    avg_z / (double) num_vertex_per_elem}.data(),
                &mElementalKdTreeData[i]);
    }
  }

}


void ExodusModel::readNodalVariables() {

  // Get variables names.
  exodusError(ex_get_var_param(
      mExodusId, "n", &mNumberNodalVariables),
              "ex_get_var_param");
  char *nm[mNumberNodalVariables];
  for (auto i = 0; i < mNumberNodalVariables; i++) { nm[i] = (char *) calloc((MAX_STR_LENGTH + 1), sizeof(char)); }

  try {
    exodusError(ex_get_variable_names(
        mExodusId, EX_NODAL, mNumberNodalVariables, nm),
                "ex_get_var_names");
    for (auto i = 0; i < mNumberNodalVariables; i++) { mNodalVariableNames.push_back(std::string(nm[i])); }

    // Get variable values.
    int time_step = 1;
    std::vector<double> buffer(mNumberVertices);
    for (auto i = 0; i < mNumberNodalVariables; i++) {
      exodusError(ex_get_nodal_var(
          mExodusId, time_step, (i + 1), mNumberVertices, buffer.data()),
                  "ex_get_nodal_var");
      mNodalVariables.insert(mNodalVariables.end(), buffer.begin(), buffer.end());
    }
  } catch (salvus_warning &e) {
    LOG() << "Reading nodal variables in file " + mExodusFileName + " raised a warning. \n"
        "This usually is not a problem, and just means the nodal variables were not defined in "
        "the file. If you don't specifically require nodal variables, you can ignore this.\n";
  }

  for (int i = 0; i < mNumberNodalVariables; i++) { free(nm[i]); }

}

void ExodusModel::readElementalVariables() {

  // Get variable names.
  exodusError(ex_get_var_param(
      mExodusId, "e", &mNumberElementalVariables), "ex_get_var_param");
  char *nm[mNumberElementalVariables];
  for (auto i = 0; i < mNumberElementalVariables; i++) {
    nm[i] = (char *) calloc((MAX_STR_LENGTH + 1), sizeof(char));
  }
  try {
    exodusError(ex_get_variable_names(mExodusId, EX_ELEM_BLOCK, mNumberElementalVariables, nm),
                "ex_get_variable_names");
    for (auto i = 0; i < mNumberElementalVariables; i++) { mElementalVariableNames.push_back(std::string(nm[i])); }

    std::vector<double> buffer(mNumberElements);
    for (auto i = 0; i < mNumberElementalVariables; i++) {
      exodusError(ex_get_var(mExodusId, 1, EX_ELEM_BLOCK, (i + 1), 1, mNumberElements, buffer.data()),
                  "ex_get_var " + mElementalVariableNames[i]);
      mElementalVariables.insert(mElementalVariables.end(), buffer.begin(), buffer.end());
    }
  } catch (salvus_warning &e) {
    LOG() << "Reading elemental variables in file " + mExodusFileName + " raised a warning. \n"
        "This usually is not a problem, and just means the elemental variables were not defined in "
        "the file. If you don't specifically require elemental variables, you can ignore this.\n";
  }

  for (int i = 0; i < mNumberElementalVariables; i++) { free (nm[i]); }

}

PetscReal ExodusModel::getNodalParameterAtNode(const std::vector<PetscReal> point,
                                               const std::string parameter_name) {

  if (!mNodalVariables.size()) {
    throw std::runtime_error(
        "You've tried to query a nodal parameter, but none are defined. "
            "Perhaps you meant to try an elemental parameter?" );
  }

  // Get spatial index.
  kdres *set = kd_nearest(mNodalKdTree, point.data());
  auto spatial_index = *(int *) kd_res_item_data(set);
  kd_res_free(set);

  // Get parameter index.
  int i = 0;
  int parameter_index;
  for (auto &name: mNodalVariableNames) {
    if (name == parameter_name) parameter_index = i;
    i++;
  }

  return mNodalVariables[parameter_index * mNumberVertices + spatial_index];

}

void ExodusModel::readConnectivity() {

  int block_id = 0;
  mElementConnectivity.resize(mVerticesPerElementPerBlock[block_id] * mNumberElements);
  exodusError(ex_get_elem_conn(mExodusId, mElementBlockIds[block_id], mElementConnectivity.data()), "ex_get_elem_conn");

}

PetscScalar ExodusModel::getElementalMaterialParameterAtVertex(const Eigen::VectorXd &elem_center,
                                                               std::string parameter_name,
                                                               const PetscInt vertex_num) const {
  assert(elem_center.size() == mNumberDimension);

  if (!mElementalVariables.size()) {
    throw std::runtime_error(
        "You've tried to query a elemental parameter, but none are defined. "
            "Perhaps you meant to try a nodal parameter?" );
  }

  // Get elemental spatial index.
  auto *set = kd_nearest(mElementalKdTree, elem_center.data());
  auto spatial_index = *(int *) kd_res_item_data(set);
  kd_res_free(set);


  /* We can change the requested parameter name in some cases. Still should eventuall change. */
  try {
    if ((std::find(mElementalVariableNames.begin(), mElementalVariableNames.end(), parameter_name + "_0") ==
        mElementalVariableNames.end())) {
      if ((parameter_name == "VP") && (std::find(
          mElementalVariableNames.begin(), mElementalVariableNames.end(), "VPV_0")
          != mElementalVariableNames.end())) {
        parameter_name = "VPV";
        throw salvus_warning("Isotropic VP requested, but can only find anisotropic VP. Using VPV!");
      } else if ((parameter_name == "VPV") && (std::find(
          mElementalVariableNames.begin(), mElementalVariableNames.end(), "VP_0")
          != mElementalVariableNames.end())) {
        parameter_name = "VP";
        throw salvus_warning("Anisotropic VPV requested, but can only find isotropic VP. Using VP!");
      } else {
        throw std::runtime_error("Requested parameter " + parameter_name + " which is not in the "
            "model file " + mExodusFileName);
      }
    }
  } catch (salvus_warning &e) {
    LOG() << e.what();
  }

  /* By the time we get here, we know that the parameter exists. */
  PetscInt parameter_index = 0; PetscInt i = 0;
  std::string full_name = parameter_name + "_" + std::to_string(vertex_num);
  for (auto &name: mElementalVariableNames) {
    if (name == full_name) { parameter_index = i; }
    i++;
  }

  return mElementalVariables[parameter_index * mNumberElements + spatial_index];
}

std::string ExodusModel::getElementType(const Eigen::VectorXd &elem_center) {

  assert(elem_center.size() == mNumberDimension);

  // Get spatial element index.
  auto *set = kd_nearest(mElementalKdTree, elem_center.data());
  auto spatial_index = *(int *) kd_res_item_data(set);
  kd_res_free(set);

  auto iterator = std::find(
      mElementalVariableNames.begin(),
      mElementalVariableNames.end(),
      "fluid");
  if (iterator == mElementalVariableNames.end()) {
    throw std::runtime_error("Fluid flag not found in mesh. Salvus currently requires a flag "
                                 "specifiying whether a given element is a fluid.");
  }

  /* Cast the float to an int :( */
  PetscInt index = iterator - mElementalVariableNames.begin();
  PetscInt type = mElementalVariables[index * mNumberElements + spatial_index];

  /* Decides what type the element is. Right now can only be elastic or fluid. */
  switch (type) {

    case 0: /* Elastic */

      switch (mNumberDimension) {

        case 2:
          return "2delastic";

        case 3:
          return "3delastic";

        default:
          throw std::runtime_error("Error on mesh dimension.");

      }

    case 1: /* Fluid */

      return "fluid";

    default:
      throw std::runtime_error("Unknown element physics. Currently we support "
                                   "[ fluid, elastic2d, elastic3d ]");

  }
}

void ExodusModel::readSideSets() {

  char *nm[mNumberSideSets];
  for (int i = 0; i < mNumberSideSets; i++) { nm[i] = (char *) calloc((MAX_STR_LENGTH + 1), sizeof(char)); }
  exodusError(ex_get_names(mExodusId, EX_SIDE_SET, nm), "ex_get_names");
  for (int i = 0; i < mNumberSideSets; i++) { mSideSetNames.push_back(std::string(nm[i])); }
  for (int i = 0; i < mNumberSideSets; i++) { free (nm[i]); }

}
std::string ExodusModel::SideSetName(const PetscInt side_set_num) {

  if (side_set_num >= mNumberSideSets) {
    throw std::runtime_error(
        "Side set " + std::to_string(side_set_num) + " is not in Exodus file. "
            "Defined side sets range from 0 to " + std::to_string(mNumberSideSets));
  }

  return mSideSetNames[side_set_num];

}

