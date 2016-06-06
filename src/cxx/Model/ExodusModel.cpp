#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Utilities/Utilities.h>
#include <stdexcept>

ExodusModel::ExodusModel(Options options) {
  mExodusFileName = options.ExodusModelFile();
  mNumberElementVertex = 4;
}

void ExodusModel::initializeParallel() {

    // Do all reading from exodus file on rank 0.
    if (MPI::COMM_WORLD.Get_rank() == 0) {
        getInitialization();
        readConnectivity();
        readCoordinates();
        readElementalVariables();
        readSideSets();
    }

  // Broadcast all scalars.
  mNumberDimension = utilities::broadcastInt(mNumberDimension);
  mNumberVertices = utilities::broadcastInt(mNumberVertices);
  mNumberElements = utilities::broadcastInt(mNumberElements);
  mNumberElementBlocks = utilities::broadcastInt(mNumberElementBlocks);
  mNumberNodeSets = utilities::broadcastInt(mNumberNodeSets);
  mNumberSideSets = utilities::broadcastInt(mNumberSideSets);

  // Broadcast all vectors.
  mNodalX = utilities::broadcastStdVecFromRoot(mNodalX);
  mNodalY = utilities::broadcastStdVecFromRoot(mNodalY);
  mVerticesPerElementPerBlock = utilities::broadcastStdVecFromRoot(mVerticesPerElementPerBlock);
  mElementConnectivity = utilities::broadcastStdVecFromRoot(mElementConnectivity);
  mElementalVariables = utilities::broadcastStdVecFromRoot(mElementalVariables);
  mElementalVariableNames = utilities::broadcastStringVecFromRank(mElementalVariableNames, 0);
  mSideSetNames = utilities::broadcastStringVecFromRank(mSideSetNames, 0);

  // Broadcast dimension specific components.
  if (mNumberDimension > 2) { mNodalZ = utilities::broadcastStdVecFromRoot(mNodalZ); }

  // Create KdTree on each processor.
  createElementalKdTree();


}

void ExodusModel::exodusError(const int retval, std::string func_name) {

  if (retval) {
    throw std::runtime_error("Error in exodus function: " + func_name + " with retval " + std::to_string((long long) retval));
  }

}

void ExodusModel::readCoordinates() {

  try {
    mNodalX.resize(mNumberVertices);
    mNodalY.resize(mNumberVertices);
    if (mNumberDimension > 2) { mNodalZ.resize(mNumberVertices); }
    if (mNumberDimension == 2) {
      exodusError(ex_get_coord(
                               mExodusId, mNodalX.data(), mNodalY.data(), NULL),
                  "ex_get_coord");
    } else {
      exodusError(ex_get_coord(
                               mExodusId, mNodalX.data(), mNodalY.data(), mNodalZ.data()),
                  "ex_get_coord");
    }
  } catch (std::exception &e) {
    utilities::print_from_root_mpi(e.what());
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }

}

void ExodusModel::getInitialization() {

  int io_ws = 0;
  int comp_ws = 8;
  try {
    mExodusId = ex_open(mExodusFileName.c_str(), EX_READ, &comp_ws, &io_ws, &mExodusVersion);
    if (mExodusId < 0) { throw std::runtime_error("Error opening exodus model file."); }

        exodusError(ex_get_init(
                            mExodusId, mExodusTitle, &mNumberDimension, &mNumberVertices, &mNumberElements,
                            &mNumberElementBlocks, &mNumberNodeSets, &mNumberSideSets),
                    "ex_get_init");
        mElementBlockIds.resize(mNumberElementBlocks);
        mVerticesPerElementPerBlock.resize(mNumberElementBlocks);
        exodusError(ex_get_elem_blk_ids(mExodusId,mElementBlockIds.data()),
                    "ex_get_elem_blk_ids");
        for(int i=0;i<mNumberElementBlocks;i++) {
          char elem_type[256];
          int num_elem_this_blk;
          int num_nodes_per_elem;
          int num_attr;
          exodusError(ex_get_elem_block	(	mExodusId,
                                                mElementBlockIds[i],
                                                elem_type,&num_elem_this_blk,
                                                &num_nodes_per_elem,
                                                &num_attr
                                                ), "ex_get_elem_block");
          mVerticesPerElementPerBlock[i] = num_nodes_per_elem;
        }
        

  } catch (std::exception &e) {
    utilities::print_from_root_mpi(e.what());
    MPI_Abort(PETSC_COMM_WORLD, -1);
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

  } else if(mNumberDimension == 3) {
      mElementalKdTree = kd_create(3);
      for (auto i = 0; i < mNumberElements; i++) {

        // Calculate element centers.
        double avg_x = 0, avg_y = 0, avg_z = 0;
        for (auto j = num_vertex_per_elem*i; j < num_vertex_per_elem*(i+1); j++) {
          avg_x += mNodalX[mElementConnectivity[j]-1];
          avg_y += mNodalY[mElementConnectivity[j]-1];
          avg_z += mNodalZ[mElementConnectivity[j]-1];
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
  exodusError(ex_get_var_names(
      mExodusId, "N", mNumberNodalVariables, nm),
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
}

PetscReal ExodusModel::getMaterialParameterAtPoint(const std::vector<double> point,
                                                   const std::string parameter_name) {

  // Ensure dimensions are consistent.
  assert(point.size() == mNumberDimension);

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
  mElementConnectivity.resize(mVerticesPerElementPerBlock[block_id]*mNumberElements);
  exodusError(ex_get_elem_conn(mExodusId, mElementBlockIds[block_id], mElementConnectivity.data()), "ex_get_elem_conn");

}

/* TODO: I am ashamed of this function. I'll change it once we standardize our paramters -mike. */
double ExodusModel::getElementalMaterialParameterAtVertex(const Eigen::VectorXd &elem_center,
                                                          std::string parameter_name,
                                                          const int vertex_num) const {
  assert(elem_center.size() == mNumberDimension);

  // Get elemental spatial index.
  auto *set = kd_nearest(mElementalKdTree, elem_center.data());
  auto spatial_index = *(int *) kd_res_item_data(set);
  kd_res_free(set);

  // Get parameter index.
  bool found = false;
  int parameter_index;
  while (!found) {
    int i = 0;
    std::string full_name = parameter_name + "_" + std::to_string((long long) vertex_num);
    for (auto &name: mElementalVariableNames) {
      if (name == full_name) {
        parameter_index = i; found = true;
      }
      i++;
    } if (!found && parameter_name == "VP") { parameter_name = "VPV"; }
  }

  return mElementalVariables[parameter_index * mNumberElements + spatial_index];
}

std::string ExodusModel::getElementType(const Eigen::VectorXd &elem_center) {

//  assert(elem_center.size() == mNumberDimension);

  // Get spatial element index.
  auto *set = kd_nearest(mElementalKdTree, elem_center.data());
  auto spatial_index = * (int *) kd_res_item_data(set);
  kd_res_free(set);

  int i = 0;
  int parameter_index = -1;
  auto full_name = "fluid";
  for (auto &name : mElementalVariableNames) {
    if (name == full_name) { parameter_index = i; }
    i++;
  }

  if(parameter_index == -1) {
    std::cout << "FLUID DEFAULT!!!!!!!" << std::endl;
    return "fluid";
    throw std::runtime_error("ERROR: `fluid` field not found in mesh!");
  }
  
  auto type = mElementalVariables[parameter_index * mNumberElements + spatial_index];
  if (type == 1) {
    return "fluid";
  } else {
    return "2delastic";
  }
}


void ExodusModel::readElementalVariables() {

  // Get variable names.
  try {
    exodusError(ex_get_var_param(
        mExodusId, "e", &mNumberElementalVariables), "ex_get_var_param");
    char *nm[mNumberElementalVariables];
    for (auto i = 0; i < mNumberElementalVariables; i++) {
      nm[i] = (char *) calloc((MAX_STR_LENGTH + 1), sizeof(char));
    }
    exodusError(ex_get_variable_names(mExodusId, EX_ELEM_BLOCK, mNumberElementalVariables, nm),
                "ex_get_variable_names");
    for (auto i = 0; i < mNumberElementalVariables; i++) { mElementalVariableNames.push_back(std::string(nm[i])); }

    int time_step = 1;
    std::vector<double> buffer(mNumberElements);
    for (auto i = 0; i < mNumberElementalVariables; i++) {
      exodusError(ex_get_var(mExodusId, 1, EX_ELEM_BLOCK, (i+1), 1, mNumberElements, buffer.data()),
          "ex_get_var " + mElementalVariableNames[i]);
      mElementalVariables.insert(mElementalVariables.end(), buffer.begin(), buffer.end());
    }
  } catch (std::runtime_error &e) {
    std::cout << e.what() << std::endl;
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }

}

void ExodusModel::readSideSets() {

  char *nm[mNumberSideSets];
  for (int i = 0; i < mNumberSideSets; i++) { nm[i] = (char *) calloc((MAX_STR_LENGTH + 1), sizeof(char)); }
  exodusError(ex_get_names(mExodusId, EX_SIDE_SET, nm), "ex_get_names");
  for (int i = 0; i < mNumberSideSets; i++) { mSideSetNames.push_back(std::string(nm[i])); }

}

