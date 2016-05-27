#include <petscviewerhdf5.h>

#include <mpi.h>
#include <fstream>
#include <Mesh/Mesh.h>
#include <Mesh/ScalarNewmark2D.h>
#include <Mesh/ElasticNewmark2D.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Utilities/Utilities.h>

void exodusError(const int retval, std::string func_name) {

  if (retval) {
    throw std::runtime_error("Error in exodus function: " + func_name);
  }

}

extern "C" {
#include <exodusII.h>
}

Mesh *Mesh::factory(Options options) {

  std::string mesh_type(options.MeshType());
  try {
    if (mesh_type == "newmark") {
      auto sc_nm_mesh = new ScalarNewmark2D();
      if (options.TestIC()) {
        // add nodal location to global field vectors for testing
        sc_nm_mesh->AddToGlobalFields("x");
        sc_nm_mesh->AddToGlobalFields("z");
        if (options.ElementShape() == "hex") {
          sc_nm_mesh->AddToGlobalFields("y");
        }
      }
      return sc_nm_mesh;
    } else if (mesh_type == "newmark_2d_elastic") {
      return new ElasticNewmark2D();
    } else {
      throw std::runtime_error("Runtime Error: Mesh type " + mesh_type + " not supported");
    }
  } catch (std::exception &e) {
    utilities::print_from_root_mpi(e.what());
    MPI_Abort(PETSC_COMM_WORLD, -1);
    return nullptr;
  };

}

#undef __FUNCT__
#define __FUNCT__ "readBoundaryNames"
int Mesh::readBoundaryNames(Options options) {

  PetscFunctionBegin;
  PetscInt rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  std::vector<std::string> boundary_names;
  if (rank == 0) {

    int num_boundaries = -1;
    float exodus_version;
    int io_ws = 0;
    int comp_ws = 8;
    try {
      int exoid = ex_open(options.ExodusMeshFile().c_str(), EX_READ, &comp_ws, &io_ws, &exodus_version);
      if (exoid < 0) { throw std::runtime_error("Error opening exodus model file."); }
      float fdum = 0.0;
      char *cdum = 0;
      exodusError(ex_inquire(exoid, EX_INQ_SIDE_SETS, &num_boundaries, &fdum, cdum), "ex_inq");
      for (int b = 1; b <= num_boundaries; b++) {
        char name[1024];
        int ret = -1;
        // retrieve name from id={1,2,3..}
        // note that blank spaces given in Trelis get an underscore.
        exodusError(ex_get_name(exoid, EX_SIDE_SET, b, name), "ex_get_name");
        std::string namestr(name);
        boundary_names.push_back(name);
      }
    } catch (std::exception &e) {
      std::cerr << "ERROR(rank=" << rank << ")!: " << e.what();
      MPI_Abort(PETSC_COMM_WORLD, -1);
    }
  } // rank==0

  boundary_names = utilities::broadcastStringVecFromRank(boundary_names, 0);

  // Build mapping boundary name -> label value (id). `id` starts counting at 1.
  for (int i = 0; i < boundary_names.size(); i++) {
    mBoundaryIds[i + 1] = boundary_names[i];
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "read"
void Mesh::read(Options options) {

  // Class variables.
  mDistributedMesh = NULL;
  mExodusFileName = options.ExodusMeshFile();

  // check if file exists
  PetscInt rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) {
    std::ifstream f(mExodusFileName.c_str());
    if (!f.good()) {
      std::cerr << "ERROR: Mesh " << mExodusFileName << " does not exist!\n";
      MPI_Abort(PETSC_COMM_WORLD, -1);
    }
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  // Function variables.
  DM dm = NULL;
  PetscBool interpolate_edges = PETSC_TRUE;

  // Read exodus file.
  DMPlexCreateExodusFromFile(PETSC_COMM_WORLD, mExodusFileName.c_str(), interpolate_edges, &dm);

  // May be a race condition on distribute.
  MPI_Barrier(PETSC_COMM_WORLD);
  DMPlexDistribute(dm, 0, NULL, &mDistributedMesh);

  // We don't need the serial mesh anymore.
  if (mDistributedMesh) { DMDestroy(&dm); }
    // mDistributedMesh == NULL when only 1 proc is used.
  else { mDistributedMesh = dm; }

  DMGetDimension(mDistributedMesh, &mNumberDimensions);
  DMPlexGetDepthStratum(mDistributedMesh, mNumberDimensions, NULL, &mNumberElementsLocal);
}

void Mesh::read(int dim, int numCells, int numVerts, int numVertsPerElem,
                Eigen::MatrixXi cells, Eigen::MatrixXd vertex_coords) {

  // Class variables.
  mDistributedMesh = NULL;
  DM dm = NULL;
  PetscBool interpolate_edges = PETSC_TRUE;

  printf("cells=[");
  for (int i = 0; i < numCells * numVertsPerElem; i++) {
    printf("%d, ", cells.data()[i]);
  }
  printf("]\n");
  DMPlexCreateFromCellList(PETSC_COMM_WORLD, dim, numCells, numVerts, numVertsPerElem,
                           interpolate_edges, cells.data(), dim, vertex_coords.data(), &dm);

  // May be a race condition on distribute.
  MPI_Barrier(PETSC_COMM_WORLD);
  DMPlexDistribute(dm, 0, NULL, &mDistributedMesh);

  // We don't need the serial mesh anymore.
  if (mDistributedMesh) { DMDestroy(&dm); }
    // mDistributedMesh == NULL when only 1 proc is used.
  else { mDistributedMesh = dm; }

  DMGetDimension(mDistributedMesh, &mNumberDimensions);
  DMPlexGetDepthStratum(mDistributedMesh, mNumberDimensions, NULL, &mNumberElementsLocal);
}

#undef __FUNCT__
#define __FUNCT__ "setupBoundaries"
int Mesh::setupBoundaries(Options options) {
  PetscFunctionBegin;

  int ierr = 0;

  readBoundaryNames(options);

  // PETSc gives all side sets the label name "Face Sets" Each face
  // is then given value, which is the order found in the exodus
  // file; Ex: two boundaries "absorbing_boundaries", "free_surface"
  // faces on absorbing boundary get label value "1", and faces on
  // free surface get value "2".
  DMLabel label;
  ierr = DMPlexGetLabel(mDistributedMesh, "Face Sets", &label);
  CHKERRQ(ierr);

  // if no side sets are present, proceed with no boundaries (naturally a free surface)
  if (label == NULL) { return 0; }
  ierr = DMLabelGetNumValues(label, &mNumberSideSets);
  CHKERRQ(ierr);
  if (mNumberSideSets == 0) { return 0; }

  IS idIS;
  const PetscInt *ids;
  ierr = DMLabelGetValueIS(label, &idIS);
  CHKERRQ(ierr);
  ierr = ISGetIndices(idIS, &ids);
  CHKERRQ(ierr);

  for (int i = 0; i < mNumberSideSets; i++) {

    IS pointIS;
    const PetscInt *faces;
    PetscInt numFaces, p;

    ierr = DMLabelGetStratumSize(label, ids[i], &numFaces);
    CHKERRQ(ierr);
    ierr = DMLabelGetStratumIS(label, ids[i], &pointIS);
    CHKERRQ(ierr);
    ierr = ISGetIndices(pointIS, &faces);
    CHKERRQ(ierr);
    auto boundary_name = mBoundaryIds[ids[i]];
    for (int p = 0; p < numFaces; p++) {

      PetscInt support_size;
      const PetscInt *support;
      int face = faces[p];
      DMPlexGetSupportSize(mDistributedMesh, face, &support_size);
      DMPlexGetSupport(mDistributedMesh, face, &support);
      for (int elem = 0; elem < support_size; elem++) {
        mBoundaryElementFaces[boundary_name][support[elem]].push_back(face);
      }
    }
    ISRestoreIndices(pointIS, &faces);
    ISDestroy(&pointIS);
  }
  ierr = ISRestoreIndices(idIS, &ids);
  CHKERRQ(ierr);
  ISDestroy(&idIS);

  PetscFunctionReturn(0);
}

PetscErrorCode Mesh::setupGlobalDof(int num_dof_vtx, int num_dof_edg,
                                    int num_dof_fac, int num_dof_vol,
                                    int num_dim, ExodusModel *model) {

  assert(num_dim == mNumberDimensions);

  // Generic error code.
  PetscErrorCode ier;

  // Resize vector of all element types.
  mElmFields.resize(mNumberElementsLocal);

  // Save the types for all the elements. Get the types with our (could be better)
  // centroid method.
  for (PetscInt i = 0; i < mNumberElementsLocal; i++) {

    // Get the type of this element.
    Eigen::MatrixXd vtx = getElementCoordinateClosure(i);
    Eigen::VectorXd ctr(vtx.cols());

    // Get the transitive closure for this element.
    PetscInt num_pnts;
    PetscInt *points = NULL;
    ier = DMPlexGetTransitiveClosure(mDistributedMesh, i, PETSC_TRUE, &num_pnts, &points);
    CHKERRQ(ier);

    // Extract the string defining the element type.
    if (mNumberDimensions == 2) {
      ctr << vtx.col(0).mean(), vtx.col(1).mean();
    }
    else if (mNumberDimensions == 3) {
      ctr << vtx.col(0).mean(), vtx.col(1).mean(), vtx.col(2).mean();
    }
    std::string type = model->getElementType(ctr);

    // Set primary fields.
    for (auto &field: appendPhysicalFields({}, type)) {
      mElmFields[i].insert(field);
      mMeshFields.insert(field);
    }

    // Traverse the closure, and add types accordingly.
    for (PetscInt j = 0; j < num_pnts; j++) {
      // Make a map of pnt -> set of types at each graph point.
      PetscInt pnt = points[2 * j];
      for (auto &field: appendPhysicalFields({}, type)) { mCouplingFields[pnt].insert(field); }
    }
  }

  // Create the section and add the fields we've defined.
  ier = PetscSectionCreate(PetscObjectComm((PetscObject) mDistributedMesh), &mMeshSection);
  CHKERRQ(ier);
  ier = PetscSectionSetNumFields(mMeshSection, mMeshFields.size());
  CHKERRQ(ier);

  int itr = 0;
  int one_field = 1;
  for (auto field: mMeshFields) {
    ier = PetscSectionSetFieldComponents(mMeshSection, itr, one_field);
    CHKERRQ(ier);
    ier = PetscSectionSetFieldName(mMeshSection, itr, field.c_str());
    CHKERRQ(ier);
    itr++;
  }

  // Get the indices of  the entire Hasse Diagram.
  PetscInt p_start, p_end;
  ier = DMPlexGetChart(mDistributedMesh, &p_start, &p_end);
  CHKERRQ(ier);

  // Tell the section about its' bounds.
  ier = PetscSectionSetChart(mMeshSection, p_start, p_end);
  CHKERRQ(ier);

  // Get the number of levels of the Hasse diagram. Should be dimension + 1 (fac, edg, vtx. for 2D)
  PetscInt depth;
  ier = DMPlexGetDepth(mDistributedMesh, &depth);
  CHKERRQ(ier);

  // Wat
  PetscInt *p_max;
  ier = PetscMalloc1(depth + 1, &p_max);
  CHKERRQ(ier);
  ier = DMPlexGetHybridBounds(mDistributedMesh,
                              depth >= 0 ? &p_max[depth] : NULL,
                              depth > 1 ? &p_max[depth - 1] : NULL,
                              depth > 2 ? &p_max[1] : NULL,
                              &p_max[0]);
  CHKERRQ(ier);

  // dep = 0 -> fac; dep = 1 -> edg; dep = 2 -> vtx.
//  assert(depth == 2); // Can only handle 2-D elements for now.
  for (int dep = 0; dep <= depth; dep++) {
    int num_dof = 0;
    if (dep == 0) { num_dof = num_dof_vtx; }
    else if (dep == 1) { num_dof = num_dof_edg; }
    else if (dep == 2) { num_dof = num_dof_fac; }
    else if (dep == 3) { num_dof = num_dof_vol; }
    PetscInt d = mNumberDimensions == depth ? dep : (!dep ? 0 : mNumberDimensions);
    // Get the number of components in each level of the Hasse diagram.
    ier = DMPlexGetDepthStratum(mDistributedMesh, dep, &p_start, &p_end);
    CHKERRQ(ier);
    p_max[dep] = p_max[dep] < 0 ? p_end : p_max[dep];
    // Loop over all mesh points on this level.
    for (int p = p_start; p < p_end; ++p) {
      PetscInt tot = 0;
      for (int f = 0; f < mMeshFields.size(); ++f) {

        // Loop through all fields defined on the mesh.
        const char *nm;
        PetscSectionGetFieldName(mMeshSection, f, &nm);

        // Is this particular field defined at this mesh point?
        bool exists = mCouplingFields[p].find((std::string(nm))) != mCouplingFields[p].end();

        // If so, add num_dof points to this field. Otherwise add 0 dofs.
        // Remember: the number of components given to any particular field
        // has already been set above.
        int num_dof_field_elem = exists ? num_dof : 0;

        // Actually set the DOFs into the section.
        ier = PetscSectionSetFieldDof(mMeshSection, p, f, num_dof);
        CHKERRQ(ier);

        // Set a custom number of dofs for each field.
        tot += num_dof;
      }
      // Total number of dofs per points is a sum of all the field dofs.
      ier = PetscSectionSetDof(mMeshSection, p, tot);
      CHKERRQ(ier);
    }
  }

  // Clean up.
  ier = PetscFree(p_max);
  CHKERRQ(ier);

  // Finalize the section, and attach it to the DM.
  ier = PetscSectionSetUp(mMeshSection);
  CHKERRQ(ier);
  DMSetDefaultSection(mDistributedMesh, mMeshSection);

  // Improves performance of closure calls (for a 1.5x application perf. boost)
  DMPlexCreateClosureIndex(mDistributedMesh, mMeshSection);

  // Go through the elements once more and assign specific couplings.
  for (PetscInt i = 0; i < mNumberElementsLocal; i++) {
    PetscInt num_pnts;
    PetscInt *points = NULL;
    ier = DMPlexGetTransitiveClosure(mDistributedMesh, i, PETSC_TRUE, &num_pnts, &points);
    CHKERRQ(ier);
    for (PetscInt j = 0; j < num_pnts; j++) {
      std::vector<std::string> v(8);
      PetscInt pnt = points[2 * j];
      std::cout << "PNT: " << pnt << std::endl;
      std::set_difference(mElmFields[i].begin(), mElmFields[i].end(),
                          mCouplingFields[pnt].begin(), mCouplingFields[pnt].end(),
                          v.begin());
//      mCouplingFields[pnt] - mElmFields[i];
//      for (auto vi: v) { std::cout << vi << std::endl; }
    }

  }


  return ier;
}

void Mesh::registerFieldVectors(const std::string &name) {

  Vec field_vector_local;
  Vec field_vector_global;
  DMCreateLocalVector(mDistributedMesh, &field_vector_local);
  DMCreateGlobalVector(mDistributedMesh, &field_vector_global);

  double zero = 0.0;
  VecSet(field_vector_local, zero);
  VecSet(field_vector_global, zero);
  PetscObjectSetName((PetscObject) field_vector_global, name.c_str());

  vec_struct registrar;
  registrar.name = name;
  registrar.loc = field_vector_local;
  registrar.glb = field_vector_global;
  mFields[name] = registrar;

}

void Mesh::checkOutField(const std::string &name) {

  assert(mFields.find(name) != mFields.end());

  // Begin the MPI broadcast global -> local.
  DMGlobalToLocalBegin(mDistributedMesh, mFields[name].glb, INSERT_VALUES,
                       mFields[name].loc);

  // End the MPI broadcast global -> local.
  DMGlobalToLocalEnd(mDistributedMesh, mFields[name].glb, INSERT_VALUES,
                     mFields[name].loc);

}

Eigen::VectorXd Mesh::getFieldOnFace(const std::string &name, const int &face_number) {

  PetscScalar *val = NULL;
  PetscInt num_nodes_face = -1;
  DMPlexVecGetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                      face_number, &num_nodes_face, &val);
  Eigen::VectorXd field(num_nodes_face);
  for (auto j = 0; j < num_nodes_face; j++) { field(j) = val[j]; }
  DMPlexVecRestoreClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                          face_number, NULL, &val);
  return field;
}

void Mesh::setFieldFromFace(const std::string &name, const int face_number, const Eigen::VectorXd &field) {

  Eigen::VectorXd val(field.size());
  // map "our" nodal ordering back onto PETSC ordering
  for (auto j = 0; j < field.size(); j++) { val(j) = field(j); }
  int ierr = DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                                 face_number, val.data(), INSERT_VALUES);
  if (ierr > 0) {
    std::cerr << "Error after DMPlexVecSetClosure in setFieldFromFace" << "\n";
  }
}

void Mesh::addFieldFromFace(const std::string &name, const int face_number, const Eigen::VectorXd &field) {

  Eigen::VectorXd val(field.size());
  // map "our" nodal ordering back onto PETSC ordering
  for (auto j = 0; j < field.size(); j++) { val(j) = field(j); }
  DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                      face_number, val.data(), ADD_VALUES);
}

Eigen::VectorXd Mesh::getFieldOnPoint(int point, std::string name) {

  int num_values;
  double *values = NULL;
  DMPlexVecGetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                      point, &num_values, &values);
  Eigen::VectorXd field(num_values);
  for (auto j = 0; j < num_values; j++) { field(j) = values[j]; }
  // printf("getting %d into %d\n",closure(47),47);
  DMPlexVecRestoreClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                          point, &num_values, &values);
  return field;
}

Eigen::VectorXd Mesh::getFieldOnElement(const std::string &name, const int &element_number,
                                        const Eigen::VectorXi &closure) {

//  DMPlexGetClosureIndices()

  PetscScalar *val = NULL;
  Eigen::VectorXd field(closure.size());
  DMPlexVecGetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                      element_number, NULL, &val);
  for (auto j = 0; j < closure.size(); j++) { field(closure(j)) = val[j]; }
  DMPlexVecRestoreClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                          element_number, NULL, &val);
  return field;

}

void Mesh::setFieldFromElement(const std::string &name, const int element_number,
                               const Eigen::VectorXi &closure, const Eigen::VectorXd &field) {

  Eigen::VectorXd val(closure.size());
  // map "our" nodal ordering back onto PETSC ordering
  for (auto j = 0; j < closure.size(); j++) { val(j) = field(closure(j)); }
  int ierr = DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                                 element_number, val.data(), INSERT_VALUES);
  if (ierr > 0) {
    std::cerr << "Error after DMPlexVecSetClosure in setFieldFromElement" << "\n";
    exit(1);
  }
}

void Mesh::addFieldFromElement(const std::string &name, const int element_number,
                               const Eigen::VectorXi &closure, const Eigen::VectorXd &field) {

  Eigen::VectorXd val(closure.size());
  // map "our" nodal ordering back onto PETSC ordering
  for (auto j = 0; j < closure.size(); j++) { val(j) = field(closure(j)); }
  int ierr = DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name].loc,
                                 element_number, val.data(), ADD_VALUES);
  if (ierr > 0) {
    std::cerr << "Error after DMPlexVecSetClosure in addFieldFromElement" << "\n";
    exit(1);
  }
}

void Mesh::assembleLocalFieldToGlobal(const std::string &name) {

  assembleLocalFieldToGlobalBegin(name);
  assembleLocalFieldToGlobalEnd(name);
}

void Mesh::assembleLocalFieldToGlobalBegin(const std::string &name) {

  // Make sure the field exists in our dictionary.
  assert(mFields.find(name) != mFields.end());

  // Begin MPI broadcast local -> global.
  DMLocalToGlobalBegin(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);
}

void Mesh::assembleLocalFieldToGlobalEnd(const std::string &name) {

  // Make sure the field exists in our dictionary.
  assert(mFields.find(name) != mFields.end());

  // Begin MPI broadcast local -> global.
  DMLocalToGlobalEnd(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);
}

void Mesh::setLocalFieldToGlobal(const std::string &name) {

  // Make sure the field exists in our dictionary.
  assert(mFields.find(name) != mFields.end());

  // Do "communication". `INSERT_VALUE` will result in no communication
  DMLocalToGlobalBegin(mDistributedMesh, mFields[name].loc, INSERT_VALUES, mFields[name].glb);
  DMLocalToGlobalEnd(mDistributedMesh, mFields[name].loc, INSERT_VALUES, mFields[name].glb);

}

void Mesh::checkInFieldBegin(const std::string &name) {

  // Make sure the field exists in our dictionary.
  assert(mFields.find(name) != mFields.end());

  // Begin MPI broadcast local -> global.
  DMLocalToGlobalBegin(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);
}


void Mesh::checkInFieldEnd(const std::string &name) {

  // Make sure the field exists in our dictionary.
  assert(mFields.find(name) != mFields.end());

  // Begin MPI broadcast local -> global.
  DMLocalToGlobalEnd(mDistributedMesh, mFields[name].loc, ADD_VALUES, mFields[name].glb);

}

void Mesh::zeroFields(const std::string &name) {
  double zero = 0.0;
  VecSet(mFields[name].loc, zero);
  VecSet(mFields[name].glb, zero);
}

void Mesh::setUpMovie(const std::string &movie_filename) {
  mViewer = nullptr;
  PetscViewerHDF5Open(PETSC_COMM_WORLD, movie_filename.c_str(), FILE_MODE_WRITE, &mViewer);
  PetscViewerHDF5PushGroup(mViewer, "/");
  DMView(mDistributedMesh, mViewer);
}

void Mesh::saveFrame(std::string name, PetscInt timestep) {

  DMSetOutputSequenceNumber(mDistributedMesh, timestep, timestep);
  int ierr = VecView(mFields[name].glb, mViewer);
  if (ierr > 0) {
    std::cerr << "ERROR @ saveFrame->VecView()\n";
    exit(1);
  }
  double max;
  VecMax(mFields[name].glb, NULL, &max);
}

void Mesh::finalizeMovie() {

  PetscViewerHDF5PopGroup(mViewer);
  PetscViewerDestroy(&mViewer);

}

Eigen::MatrixXd Mesh::getElementCoordinateClosure(PetscInt elem_num) {

  Vec coord;
  DMGetCoordinatesLocal(mDistributedMesh, &coord);
  PetscSection coord_section;
  DMGetCoordinateSection(mDistributedMesh, &coord_section);

  PetscInt coord_buf_size;
  PetscReal *coord_buf = NULL;
  DMPlexVecGetClosure(mDistributedMesh, coord_section, coord, elem_num, &coord_buf_size, &coord_buf);

  int num_vtx = coord_buf_size / mNumberDimensions;
  Eigen::MatrixXd vtx(num_vtx, mNumberDimensions);
  for (int i = 0; i < num_vtx; i++) {
    vtx(i, 0) = coord_buf[mNumberDimensions * i + 0];
    vtx(i, 1) = coord_buf[mNumberDimensions * i + 1];
    if (mNumberDimensions == 3) {
      vtx(i, 2) = coord_buf[mNumberDimensions * i + 2];
    }
  }

  DMPlexVecRestoreClosure(mDistributedMesh, coord_section, coord, elem_num, &coord_buf_size, &coord_buf);

  return vtx;

}

std::vector<std::string> Mesh::appendPhysicalFields(const std::vector<std::string>& fields,
                                                    const std::string& physics) {

  std::vector<std::string> new_fields;
  if (physics == "fluid") {
    new_fields = {"u"};
  } else if (physics == "2delastic") {
    new_fields = {"ux", "uy"};
  } else {
    PRINT_ROOT() << "ERROR IN FIELD TYPE " + physics;
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }

  for (auto &f: fields) { new_fields.push_back(f); }
  return new_fields;

}



int Mesh::numFieldPerPhysics(std::string physics) {
  int num;
  try {
    if (physics == "fluid") { num = 1; }
    else if (physics == "2delastic") { num = 2; }
    else {
      throw std::runtime_error("Derived type " + physics + " is not known.");
    }
  } catch (std::exception &e) {
    PRINT_ROOT() << e.what();
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }
  return num;
}

// int Mesh::BoundarylementFaces(int elm, int ss_num) {

//     // Make sure we're looking for a side set that actually exists.
//     // assert(mBoundaryElementFaces.find(ss_num) != mBoundaryElementFaces.end());

//     // if (mBoundaryElementFaces[ss_num].find(elm) == mBoundaryElementFaces[ss_num].end()) {
//     //     return -1;
//     // } else {
//     //     return mBoundaryElementFaces[ss_num][elm][0];
//     // }

// }

