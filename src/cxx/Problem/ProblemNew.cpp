#include <Problem/ProblemNew.h>
#include <Utilities/Options.h>
#include <Utilities/Logging.h>
#include <Problem/Order2Newmark.h>
#include <Mesh/Mesh.h>
#include <stdexcept>

using namespace Eigen;

std::unique_ptr<ProblemNew> ProblemNew::Factory(std::unique_ptr<Options> const &options) {

  std::unique_ptr<ProblemNew> problem;
  std::string timestep_scheme = "newmark";
  try {

    if (timestep_scheme == "newmark") {
      return std::unique_ptr<ProblemNew> (new Order2Newmark);
    }
    else
    {
      throw std::runtime_error("Runtime Error. Problem type not defined.");
    }

  } catch (std::runtime_error e) {
    MPI_Abort(PETSC_COMM_WORLD, -1);
    LOG() << e.what();
  }

  return problem;

}

ElemVec ProblemNew::initializeElements(unique_ptr<Mesh> const &mesh,
                                       unique_ptr<ExodusModel> const &model,
                                       unique_ptr<Options> const &options) {

  /* MPI rank is important for source/receiver detection. */
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  /* Define variables to test for source/receiver presence. */
  bool true_attach = true;
  bool trial_attach = false;

  /* Create a vector of sources and receivers from user options. */
  auto srcs = Source::Factory(options);
  auto recs = Receiver::Factory(options);

  /* Keep local track of srcs/recs on this partition. */
  std::vector<PetscInt> srcs_this_partition(srcs.size(), false);
  std::vector<PetscInt> recs_this_partition(recs.size(), false);

  /* All of our (polymorphic) elements will lie here. */
  ElemVec elements;

  /* Allocate all elements. */
  for (PetscInt i = 0; i < mesh->NumberElementsLocal(); i++)
  {

    /* Push back an appropriate element based on the mesh. */
    elements.push_back(Element::Factory(mesh->ElementFields(i),
                                        mesh->TotalCouplingFields(i),
                                        options));

    /* Assign a (processor-specific) number to this element. */
    elements.back()->SetNum(i);

    /* Attach vertex co-ordinates from mesh. */
    elements.back()->attachVertexCoordinates(mesh);

    /* Set any boundary conditions on this element. */
    elements.back()->setBoundaryConditions(mesh);

    /* Attach material properties (velocity, Cij, etc...). */
    elements.back()->attachMaterialProperties(model);

    /* Test for any sources. */
    for (auto &src: srcs) {
      srcs_this_partition[src->Num()] = elements.back()->attachSource(src, trial_attach) ? rank : 0;
    }

    /* Test for any receivers. */
    for (auto &rec: recs) {
      recs_this_partition[rec->Num()] = elements.back()->attachReceiver(rec, trial_attach) ? rank : 0;
    }

  }

  /* Check sources and receivers across all parallel partitions. */
  MPI_Allreduce(MPI_IN_PLACE, srcs_this_partition.data(), srcs_this_partition.size(),
                MPIU_INT, MPI_MAX, PETSC_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, recs_this_partition.data(), recs_this_partition.size(),
                MPIU_INT, MPI_MAX, PETSC_COMM_WORLD);

  /* Finish up adding parallel-aware sources and receivers. */
  for (auto &elm: elements) {
    for (auto &src: srcs) {
      if (!src) { continue; }
      if (srcs_this_partition[src->Num()] == rank) { elm->attachSource(src, true_attach); }
    }
    for (auto &rec: recs) {
      if (!rec) { continue; }
      if (recs_this_partition[rec->Num()] == rank) { elm->attachReceiver(rec, true_attach); }
    }
  }

  /* Finally, go back and ensure that everything has been added as expected. */
  try {
    for (auto &src: srcs) {
      /* Was there a source that should have been added by this processor that wasn't? */
      if (src && srcs_this_partition[src->Num()] == rank) {
        throw std::runtime_error("Error. One or more sources were not added properly.");
      }
    }
    for (auto &rec: recs) {
      /* Was there a receiver that should have been added by this processor that wasn't? */
      if (rec && recs_this_partition[rec->Num()] == rank) {
        throw std::runtime_error("Error. One or more receivers was not added properly.");
      }
    }
  } catch (std::exception &e) {
    LOG() << e.what();
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }

  /* Now all elements are complete. */
  return elements;

}
std::tuple<ElemVec, FieldDict> ProblemNew::assembleIntoGlobalDof(ElemVec elements, FieldDict fields,
                                                                 DM PETScDM, PetscSection PETScSection,
                                                                 std::unique_ptr<Options> const &options) {

  /* Get some derived quantities. */
  PetscInt maxLocalFields = 0;
  PetscInt numDofPerElm = elements[0]->NumIntPnt();
  std::set<std::string> pullFields, pushFields;
  for (auto &elm: elements) {
    PetscInt numLocalFields = elm->PullElementalFields().size();
    maxLocalFields = numLocalFields > maxLocalFields ? numLocalFields : maxLocalFields;
    for (auto &f: elm->PullElementalFields()) {
      pullFields.insert(f);
      pushFields.insert(f);
    }
  }

  /* Matrices to hold pulled values on each element. */
  Matrix<PetscReal, Dynamic, Dynamic> a(numDofPerElm, maxLocalFields); /// acceleration.
  Matrix<PetscReal, Dynamic, Dynamic> s(numDofPerElm, maxLocalFields); /// surf/edg integral.
  Matrix<PetscReal, Dynamic, Dynamic> f(numDofPerElm, maxLocalFields); /// forcing.
  Matrix<PetscReal, Dynamic, Dynamic> u(numDofPerElm, maxLocalFields); /// vol/face integral.
  Matrix<PetscReal, Dynamic, Dynamic> k(numDofPerElm, maxLocalFields); /// stiffness.

  /* Get fields on local partitions. */
  for (auto &field: pullFields) { checkOutField(field, PETScDM, fields); }

  /* Zero fields to which we will assemble. */
  for (auto &field: pushFields) { zeroField(field, fields); }

  /* Loop over elements and compute integrals. */
  for (auto &elm: elements) {

    /* Specific number of fields for this element. */
    PetscInt NumPullFields = elm->PullElementalFields().size();
    PetscInt NumPushFields = elm->PushElementalFields().size();

    /* Get specific fields for this element. */
    for (PetscInt i = 0; i < NumPullFields; i++) {
      u.col(i) =
          getFieldOnElement(elm->PullElementalFields()[i], elm->Num(), elm->ClsMap(),
                            PETScDM, PETScSection, fields);
    }

    /* Compute stiffness. */
    k.leftCols(NumPushFields) = elm->computeStiffnessTerm(u.leftCols(NumPullFields));

    /* Compute surface integrals. */
    s.leftCols(NumPushFields) = elm->computeSurfaceIntegral(u.leftCols(NumPullFields));

    /* Compute forcing. */
    PetscReal time = 0;
    f.leftCols(NumPushFields) = elm->computeSourceTerm(time);

    /* Compute acceleration. */
    a.leftCols(NumPushFields) = f.leftCols(NumPushFields).array() -
        k.leftCols(NumPushFields).array() + s.leftCols(NumPushFields).array();

    /* Assemble fields into local partition. */
    for (PetscInt i = 0; i < NumPushFields; i++) {
      addFieldOnElement(elm->PushElementalFields()[i], elm->Num(), elm->ClsMap(),
                        a.col(i), PETScDM, PETScSection, fields);
    }

  }

  /* Broadcast to global partitions. */
  for (auto &field: pushFields) { checkInField(field, PETScDM, fields); }

  return std::tuple<ElemVec, FieldDict> (std::move(elements), std::move(fields));

}

void ProblemNew::zeroField(const std::string &name, FieldDict &fields) {

  /* Set both local and global vectors to zero. */
  PetscReal zero = 0.0;
  VecSet(fields[name]->mLoc, zero);
  VecSet(fields[name]->mGlb, zero);

}

void ProblemNew::checkInField(const std::string &name, DM PETScDM, FieldDict &fields) {

  /* Ensure field is registered. */
  assert(fields.find(name) != fields.end());

  /* Assemble local -> global. */
  DMLocalToGlobalBegin(PETScDM, fields[name]->mLoc, ADD_VALUES, fields[name]->mGlb);
  DMLocalToGlobalEnd(PETScDM, fields[name]->mLoc, ADD_VALUES, fields[name]->mGlb);

}


void ProblemNew::checkOutField(const std::string &name, DM PETScDM, FieldDict &fields) {

  /* Ensure field is registered. */
  assert(fields.find(name) != fields.end());

  /* Broadcast global -> local. */
  DMGlobalToLocalBegin(PETScDM, fields[name]->mGlb, INSERT_VALUES, fields[name]->mLoc);
  DMGlobalToLocalEnd(PETScDM, fields[name]->mGlb, INSERT_VALUES, fields[name]->mLoc);

}

void ProblemNew::addFieldOnElement(const std::string &name,
                                   const PetscInt num,
                                   const Eigen::Ref<const IntVec> &closure,
                                   const Eigen::Ref<const RealVec> &field,
                                   DM PETScDM,
                                   PetscSection PETScSection,
                                   FieldDict &fields) {

  RealVec reorder(closure.size());
  for (PetscInt i = 0; i < closure.size(); i++) { reorder(i) = field(closure(i)); }
  DMPlexVecSetClosure(PETScDM, PETScSection, fields[name]->mLoc, num, reorder.data(), ADD_VALUES);

}

RealVec ProblemNew::getFieldOnElement(const std::string &name,
                                      const PetscInt num,
                                      const Eigen::Ref<const IntVec> &closure,
                                      DM PETScDM,
                                      PetscSection PETScSection,
                                      FieldDict &fields) {

  /* Initialize arrays to get global DOF values. RVO should apply. */
  PetscScalar *val = NULL;
  RealVec field(closure.size());

  /* Populate 'val' with field on element, in PETSc ordering. */
  DMPlexVecGetClosure(PETScDM, PETScSection, fields[name]->mLoc, num, NULL, &val);

  /* Re-order into our tensorized basis ordering. */
  for (PetscInt i = 0; i < closure.size(); i++) { field(closure(i)) = val[i]; }

  /* De-allocate the PETSc pointer. */
  DMPlexVecRestoreClosure(PETScDM, PETScSection, fields[name]->mLoc, num, NULL, &val);

  return field;

}







