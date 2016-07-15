#include <petscviewerhdf5.h>
#include <stdexcept>

#include <mpi.h>
#include <fstream>
#include <Mesh/Mesh.h>
#include <Mesh/ScalarNewmark2D.h>
#include <Mesh/ElasticNewmark2D.h>
#include <Mesh/ElasticAcousticNewmark2D.h>
#include <Mesh/ElasticAcousticNewmark3D.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Utilities/Utilities.h>
#include <Utilities/Logging.h>
#include <petscdmplex.h>
#include <Utilities/Types.h>

extern "C" {
#include <exodusII.h>
}

using std::unique_ptr;

void exodusError(const int retval, std::string func_name) {

  if (retval) {
    throw std::runtime_error("Error in exodus function: " + func_name);
  }

}

Mesh::Mesh(const std::unique_ptr<Options> &options) {
  mDistributedMesh = NULL;
}

std::unique_ptr<Mesh> Mesh::factory(const std::unique_ptr<Options> &options) {

  std::string mesh_type(options->MeshType());
  try {
    if (mesh_type == "newmark") {
      std::unique_ptr<Mesh> sc_nm_mesh(new ScalarNewmark2D(options));
      if (options->TestIC()) {
        // add nodal location to global field vectors for testing
        sc_nm_mesh->AddToGlobalFields("x");
        sc_nm_mesh->AddToGlobalFields("z");
        if (options->ElementShape() == "hex") {
          sc_nm_mesh->AddToGlobalFields("y");
        }
      }
      return sc_nm_mesh;
    } else if (mesh_type == "newmark_2d_elastic") {
      return std::unique_ptr<Mesh>(new ElasticNewmark2D(options));
    } else if (mesh_type == "2d_couple") {
      return std::unique_ptr<Mesh>(new ElasticAcousticNewmark2D(options));
    } else if (mesh_type == "3d_couple") {
      return std::unique_ptr<Mesh>(new ElasticAcousticNewmark3D(options));
    } else {
      throw std::runtime_error("Runtime Error: Mesh type " + mesh_type + " not supported");
    }
  } catch (std::exception &e) {
    LOG() << e.what();
    MPI_Abort(PETSC_COMM_WORLD, -1);
    return nullptr;
  };

}

#undef __FUNCT__
#define __FUNCT__ "readBoundaryNames"
int Mesh::readBoundaryNames(std::unique_ptr<Options> const &options) {

  PetscFunctionBegin;
  PetscInt rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  std::vector<std::string> boundary_names;
  if (rank == 0) {

    int num_boundaries = -1;
    float exodus_version;
    int io_ws = 0;
    int comp_ws = 8;
    try {
      int exoid = ex_open(options->ExodusMeshFile().c_str(), EX_READ, &comp_ws, &io_ws, &exodus_version);
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
      LOG() << "ERROR(rank=" << rank << ")!: " << e.what();
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
void Mesh::read(std::unique_ptr<Options> const &options) {

  // Class variables.
  mDistributedMesh = NULL;
  mExodusFileName = options->ExodusMeshFile();

  // check if file exists
  PetscInt rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) {
    std::ifstream f(mExodusFileName.c_str());
    if (!f.good()) {
      LOG() << "ERROR: Mesh " << mExodusFileName << " does not exist!";
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

  DMGetDimension(mDistributedMesh, &mNumDim);
  DMPlexGetDepthStratum(mDistributedMesh, mNumDim, NULL, &mNumberElementsLocal);
}

#undef __FUNCT__
#define __FUNCT__ "setupBoundaries"
int Mesh::setupBoundaries(std::unique_ptr<Options> const &options) {
  PetscFunctionBegin;

  int ierr = 0;

  readBoundaryNames(options);

  // PETSc gives all side sets the label name "Face Sets" Each face
  // is then given value, which is the order found in the exodus
  // file; Ex: two boundaries "absorbing_boundaries", "free_surface"
  // faces on absorbing boundary get label value "1", and faces on
  // free surface get value "2".
  DMLabel label;
  ierr = DMGetLabel(mDistributedMesh, "Face Sets", &label);
  CHKERRQ(ierr);

  // if no side sets are present, proceed with no boundaries (naturally a free surface)
  if (label == NULL) { return 0; }
  ierr = DMLabelGetNumValues(label, &mNumberSideSets); CHKERRQ(ierr);
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
                                    int num_dim, unique_ptr<ExodusModel> const &model) {


  /* Allocate vector which will hold element types. */
  mElmFields.resize(mNumberElementsLocal);

  /* Find all the mesh boundaries. */
  DMLabel label; DMGetLabel(mDistributedMesh, "Face Sets", &label);
  PetscInt size; DMLabelGetNumValues(label, &size);
  IS idIS; DMLabelGetValueIS(label, &idIS);
  const PetscInt *ids; ISGetIndices(idIS, &ids);
  for (PetscInt i = 0; i < size; i++) {
    PetscInt numFaces; DMLabelGetStratumSize(label, ids[i], &numFaces);
    IS pointIs; DMLabelGetStratumIS(label, ids[i], &pointIs);
    const PetscInt *faces; ISGetIndices(pointIs, &faces);
    for (PetscInt j = 0; j < numFaces; j++) { mBndPts.insert(faces[j]); }
  }

  /* Walk through the mesh and extract element types. */
  for (PetscInt i = 0; i < mNumberElementsLocal; i++) {

    /* Get and save the type of element i. */
    RealMat vtx = getElementCoordinateClosure(i); RealVec ctr(vtx.cols());
    for (PetscInt j = 0; j < mNumDim; j++) { ctr(j) = vtx.col(j).mean(); }
    std::string type = model->getElementType(ctr); mPointFields[i].insert(type);

    /* Add the type of element i to all mesh points connected via the Hasse graph. */
    PetscInt num_pts; const PetscInt *pts = NULL;
    DMPlexGetCone(mDistributedMesh, i, &pts); DMPlexGetConeSize(mDistributedMesh, i, &num_pts);
    for (PetscInt j = 0; j < num_pts; j++) {
      mPointFields[pts[j]].insert(type);
      if (mBndPts.find(pts[j]) != mBndPts.end()) {
        mPointFields[pts[j]].insert("boundary");
      }
    }

    /* Finally, add the type type to the global fields. */
    mMeshFields.insert(type);

  }

  /* Use the information extracted above to inform the global DOF layout. */
  PetscInt num_bc = 0;
  PetscInt poly_order = 4;
  PetscInt num_fields = mMeshFields.size();
  PetscInt *num_comps; PetscMalloc1(num_fields, &num_comps);
  {
    PetscInt i = 0;
    for (auto &f: mMeshFields) { num_comps[i++] = numFieldPerPhysics(f); }
  }
  PetscInt *num_dof; PetscMalloc1(num_fields * (mNumDim + 1), &num_dof);

  /* For each field... */
  for (PetscInt f = 0; f < num_fields; f++) {

    /* For each element point... (vertex, edge, face, volume) */
    for (PetscInt d = 0; d < (mNumDim + 1); d++) {

      /* Number of dofs per this element point. */
      num_dof[f * (mNumDim + 1) + d] = PetscPowInt(poly_order - 1, d) * num_comps[f];

    }
  }

  /* Generate the section which will define the global dofs. */
  DMPlexCreateSection(mDistributedMesh, mNumDim, num_fields, num_comps, num_dof,
                      num_bc, NULL, NULL, NULL, NULL, &mMeshSection);
  PetscFree(num_dof); PetscFree(num_comps);

  /* Attach some meta-information to the section. */
  {
    PetscInt i = 0;
    for (auto &f: mMeshFields) { PetscSectionSetFieldName(mMeshSection, i++, f.c_str()); }
  }

  /* Attach the section to our DM, and set the spectral ordering. */
  DMSetDefaultSection(mDistributedMesh, mMeshSection);
  DMPlexCreateSpectralClosurePermutation(mDistributedMesh, NULL);

}

std::vector<PetscInt> Mesh::EdgeNumbers(const PetscInt elm) {

  std::vector<PetscInt> edges;
  PetscInt num_pts, e_start, e_end, *points = NULL, dep_edge = 1;
  DMPlexGetTransitiveClosure(mDistributedMesh, elm, PETSC_TRUE, &num_pts, &points);
  // get range of values which define edges.
  DMPlexGetDepthStratum(mDistributedMesh, dep_edge, &e_start, &e_end);
  for (PetscInt i = 0; i < num_pts; i++) {
    PetscInt pnt = points[2*i];
    // if edge is in closure, push it back.
    if (pnt >= e_start && pnt < e_end) {
      edges.push_back(pnt);
    }
  }
  DMPlexRestoreTransitiveClosure(mDistributedMesh, elm, PETSC_TRUE, &num_pts, &points);
  return edges;
}

PetscInt Mesh::GetNeighbouringElement(const PetscInt interface, const PetscInt this_elm) const {

  PetscInt neighbour = -1, num_pts, e_start, e_end, *points = NULL;
  DMPlexGetTransitiveClosure(mDistributedMesh, interface, PETSC_FALSE, &num_pts, &points);
  for (PetscInt i = 0; i < num_pts; i++) {
    PetscInt pnt = points[2*i];
    if (pnt != interface && pnt != this_elm) {
      neighbour = pnt;
    }
  }
  DMPlexRestoreTransitiveClosure(mDistributedMesh, interface, PETSC_FALSE, &num_pts, &points);

  if (neighbour == -1) {
    LOG() << "ERROR IN GetNeighbouringElement";
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }
  return neighbour;
}

std::vector<std::tuple<PetscInt,std::vector<std::string>>> Mesh::CouplingFields(const PetscInt elm) {

  std::vector<std::tuple<PetscInt,std::vector<std::string>>> couple;
  std::set<std::string> elm_fields = mElmFields[elm];

  PetscInt num_pts, e_start, e_end;
  PetscInt *points = NULL;
  PetscErrorCode ier = DMPlexGetTransitiveClosure(mDistributedMesh, elm, PETSC_TRUE, &num_pts, &points);
  PetscInt dep_edge = 1;
  DMPlexGetDepthStratum(mDistributedMesh, dep_edge, &e_start, &e_end);
  for (PetscInt j = 0; j < num_pts; j++) {
    std::vector<std::string> field;
    PetscInt pnt = points[2*j];
    for (auto f: mPointFields[pnt]) {
      bool coupling_edge = mElmFields[elm].find(f) == mElmFields[elm].end();
      // if we're an edge (not vertex)
      if (coupling_edge && pnt >= e_start && pnt < e_end) { field.push_back(f); }
    }
    if (field.size()) {
      couple.push_back(std::make_tuple(pnt, field));
    }
  }
  DMPlexRestoreTransitiveClosure(mDistributedMesh, elm, PETSC_TRUE, &num_pts, &points);
  return couple;
}

std::vector<std::string> Mesh::TotalCouplingFields(const PetscInt elm) {

  std::set<std::string> coupling_fields;
  PetscInt num_pts; const PetscInt *pts = NULL;
  DMPlexGetCone(mDistributedMesh, elm, &pts); DMPlexGetConeSize(mDistributedMesh, elm, &num_pts);
  for (PetscInt j = 0; j < num_pts; j++) {
    for (auto &f: mPointFields[pts[j]]) {
      if (mPointFields[elm].find(f) == mPointFields[elm].end()) {
        coupling_fields.insert(f);
      }
    }
  }
  return std::vector<std::string> (coupling_fields.begin(), coupling_fields.end());

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

  unique_ptr<vec_struct> registrar(new vec_struct);
  registrar->name = name;
  registrar->loc = field_vector_local;
  registrar->glb = field_vector_global;
  mFields.insert(std::pair<std::string,unique_ptr<vec_struct>>(name, std::move(registrar)));

}

void Mesh::checkOutField(const std::string &name) {

  assert(mFields.find(name) != mFields.end());

  // Begin the MPI broadcast global -> local.
  DMGlobalToLocalBegin(mDistributedMesh, mFields[name]->glb, INSERT_VALUES,
                       mFields[name]->loc);

  // End the MPI broadcast global -> local.
  DMGlobalToLocalEnd(mDistributedMesh, mFields[name]->glb, INSERT_VALUES,
                     mFields[name]->loc);

}

Eigen::VectorXd Mesh::getFieldOnFace(const std::string &name, const int &face_number) {

  PetscScalar *val = NULL;
  PetscInt num_nodes_face = -1;
  DMPlexVecGetClosure(mDistributedMesh, mMeshSection, mFields[name]->loc,
                      face_number, &num_nodes_face, &val);
  Eigen::VectorXd field(num_nodes_face);
  for (auto j = 0; j < num_nodes_face; j++) { field(j) = val[j]; }
  DMPlexVecRestoreClosure(mDistributedMesh, mMeshSection, mFields[name]->loc,
                          face_number, NULL, &val);
  return field;
}

void Mesh::setFieldFromFace(const std::string &name, const int face_number, const Eigen::VectorXd &field) {

  Eigen::VectorXd val(field.size());
  // map "our" nodal ordering back onto PETSC ordering
  for (auto j = 0; j < field.size(); j++) { val(j) = field(j); }
  int ierr = DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name]->loc,
                                 face_number, val.data(), INSERT_VALUES);
  if (ierr > 0) {
    LOG() << "Error after DMPlexVecSetClosure in setFieldFromFace";
  }
}

void Mesh::addFieldFromFace(const std::string &name, const int face_number, const Eigen::VectorXd &field) {

  Eigen::VectorXd val(field.size());
  // map "our" nodal ordering back onto PETSC ordering
  for (auto j = 0; j < field.size(); j++) { val(j) = field(j); }
  DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name]->loc,
                      face_number, val.data(), ADD_VALUES);
}

Eigen::VectorXd Mesh::getFieldOnPoint(int point, std::string name) {

  int num_values;
  double *values = NULL;
  DMPlexVecGetClosure(mDistributedMesh, mMeshSection, mFields[name]->loc,
                      point, &num_values, &values);
  Eigen::VectorXd field(num_values);
  for (auto j = 0; j < num_values; j++) { field(j) = values[j]; }
  DMPlexVecRestoreClosure(mDistributedMesh, mMeshSection, mFields[name]->loc,
                          point, &num_values, &values);
  return field;
}

Eigen::VectorXd Mesh::getFieldOnElement(const std::string &name, const int &element_number,
                                        const Eigen::VectorXi &closure) {

  PetscScalar *val = NULL;
  Eigen::VectorXd field(closure.size());
  DMPlexVecGetClosure(mDistributedMesh, mMeshSection, mFields[name]->loc,
                      element_number, NULL, &val);
  for (auto j = 0; j < closure.size(); j++) { field(closure(j)) = val[j]; }
  DMPlexVecRestoreClosure(mDistributedMesh, mMeshSection, mFields[name]->loc,
                          element_number, NULL, &val);
  return field;

}

void Mesh::setFieldFromElement(const std::string &name, const int element_number,
                               const Eigen::VectorXi &closure, const Eigen::VectorXd &field) {

  Eigen::VectorXd val(closure.size());
  // map "our" nodal ordering back onto PETSC ordering
  for (auto j = 0; j < closure.size(); j++) { val(j) = field(closure(j)); }
  int ierr = DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name]->loc,
                                 element_number, val.data(), INSERT_VALUES);
  if (ierr > 0) {
    LOG() << "Error after DMPlexVecSetClosure in setFieldFromElement";
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }
}

void Mesh::addFieldFromElement(const std::string &name, const int element_number,
                               const Eigen::VectorXi &closure, const Eigen::VectorXd &field) {

  Eigen::VectorXd val(closure.size());
  // map "our" nodal ordering back onto PETSC ordering
  for (auto j = 0; j < closure.size(); j++) { val(j) = field(closure(j)); }
  int ierr = DMPlexVecSetClosure(mDistributedMesh, mMeshSection, mFields[name]->loc,
                                 element_number, val.data(), ADD_VALUES);
  if (ierr > 0) {
    LOG() << "Error after DMPlexVecSetClosure in addFieldFromElement";
    MPI_Abort(PETSC_COMM_WORLD, -1);
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
  DMLocalToGlobalBegin(mDistributedMesh, mFields[name]->loc, ADD_VALUES, mFields[name]->glb);
}

void Mesh::assembleLocalFieldToGlobalEnd(const std::string &name) {

  // Make sure the field exists in our dictionary.
  assert(mFields.find(name) != mFields.end());

  // Begin MPI broadcast local -> global.
  DMLocalToGlobalEnd(mDistributedMesh, mFields[name]->loc, ADD_VALUES, mFields[name]->glb);
}

void Mesh::setLocalFieldToGlobal(const std::string &name) {

  // Make sure the field exists in our dictionary.
  assert(mFields.find(name) != mFields.end());

  // Do "communication". `INSERT_VALUE` will result in no communication
  DMLocalToGlobalBegin(mDistributedMesh, mFields[name]->loc, INSERT_VALUES, mFields[name]->glb);
  DMLocalToGlobalEnd(mDistributedMesh, mFields[name]->loc, INSERT_VALUES, mFields[name]->glb);

}

void Mesh::checkInFieldBegin(const std::string &name) {

  // Make sure the field exists in our dictionary.
  assert(mFields.find(name) != mFields.end());

  // Begin MPI broadcast local -> global.
  DMLocalToGlobalBegin(mDistributedMesh, mFields[name]->loc, ADD_VALUES, mFields[name]->glb);
}


void Mesh::checkInFieldEnd(const std::string &name) {

  // Make sure the field exists in our dictionary.
  assert(mFields.find(name) != mFields.end());

  // Begin MPI broadcast local -> global.
  DMLocalToGlobalEnd(mDistributedMesh, mFields[name]->loc, ADD_VALUES, mFields[name]->glb);

}

void Mesh::zeroFields(const std::string &name) {
  double zero = 0.0;
  VecSet(mFields[name]->loc, zero);
  VecSet(mFields[name]->glb, zero);
}

void Mesh::setUpMovie(const std::string &movie_filename) {
  mViewer = nullptr;
  PetscViewerHDF5Open(PETSC_COMM_WORLD, movie_filename.c_str(), FILE_MODE_WRITE, &mViewer);
  PetscViewerHDF5PushGroup(mViewer, "/");
  DMView(mDistributedMesh, mViewer);
  int_tstep = 0;
}

void Mesh::saveFrame(std::string name, PetscInt timestep) {

  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) LOG() << "SAVING FRAME: " << name << ' ' << timestep;
  DMSetOutputSequenceNumber(mDistributedMesh, timestep, timestep);
  int ierr = VecView(mFields[name]->glb, mViewer);
  if (ierr > 0) {
    LOG() << "ERROR @ saveFrame->VecView()";
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }
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

  int num_vtx = coord_buf_size / mNumDim;
  Eigen::MatrixXd vtx(num_vtx, mNumDim);
  for (int i = 0; i < num_vtx; i++) {
    vtx(i, 0) = coord_buf[mNumDim * i + 0];
    vtx(i, 1) = coord_buf[mNumDim * i + 1];
    if (mNumDim == 3) {
      vtx(i, 2) = coord_buf[mNumDim * i + 2];
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
  } else if (physics == "3delastic") {
    new_fields = {"ux", "uy", "uz"};
  } else {
    LOG() << "ERROR IN FIELD TYPE " + physics;
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
    else if (physics == "3delastic") { num = 3; }
    else {
      throw std::runtime_error("Derived type " + physics + " is not known.");
    }
  } catch (std::exception &e) {
    LOG() << e.what();
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }
  return num;
}

