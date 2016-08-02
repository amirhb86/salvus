#include <Mesh/Mesh.h>
#include <Problem/Order2Newmark.h>
#include <Utilities/Logging.h>
#include <Utilities/Options.h>


FieldDict Order2Newmark::initializeGlobalDofs(ElemVec const &elements,
                                              std::unique_ptr<Mesh> &mesh) {

  FieldDict fields;

  /* Initialize vector which will hold diagonal mass matrix. */
  fields.insert(std::pair<std::string, std::unique_ptr<field>>
                    ("mi", std::unique_ptr<field> (new field("mi", mesh->DistributedMesh()))));

  /* Sum mass matrix into local partition. */
  for (auto &elm: elements) {
    DMPlexVecSetClosure(mesh->DistributedMesh(), mesh->MeshSection(), fields["mi"]->mLoc,
                        elm->Num(), elm->assembleElementMassMatrix().data(), ADD_VALUES);
  }

  /* Sum mass matrix into global partition. */
  DMLocalToGlobalBegin(mesh->DistributedMesh(), fields["mi"]->mLoc, ADD_VALUES, fields["mi"]->mGlb);
  DMLocalToGlobalEnd(mesh->DistributedMesh(), fields["mi"]->mLoc, ADD_VALUES, fields["mi"]->mGlb);

  /* Take component wise inverse of mass "matrix". */
  VecReciprocal(fields["mi"]->mLoc); VecReciprocal(fields["mi"]->mGlb);

  /* Initialize global field vectors. */
  if (!mesh->GlobalFields().empty()) {
    for (auto &f: mesh->GlobalFields()) {
    fields.insert(std::pair<std::string,std::unique_ptr<field>>
                      (f, std::unique_ptr<field> (new field(f, mesh->DistributedMesh()))));
    }
  } else {
    LOG() << "No global fields defined!"; MPI_Abort(PETSC_COMM_WORLD, -1);
  }

  return fields;

}

FieldDict Order2Newmark::applyInverseMassMatrix(FieldDict fields) {

  /* Ensure there is an inverted diagonal mass matrix. */
  assert(fields.find("mi") != fields.end());

  /* Multiply acceleration by inverse mass matrix (if a component exists) */
  std::vector<std::string> recognized_acl {"ax", "ay", "ax", "a"};
  for (auto &f: {"ax", "ay", "az", "a"}) {
    if (fields.count(f)) {
      VecPointwiseMult(fields[f]->mGlb, fields["mi"]->mGlb, fields[f]->mGlb);
    } else { continue; }
  }

  return fields;

}

std::tuple<FieldDict, PetscScalar> Order2Newmark::takeTimeStep(
    FieldDict fields, PetscScalar time, std::unique_ptr<Options> const &options) {

  /* Constants from Newmark scheme. */
  PetscReal acl_factor = (1.0/2.0) * mDt;
  PetscReal dsp_factor = (1.0/2.0) * (mDt * mDt);

  /* List of vectors to multiply. */
  const std::vector<std::string> recognized_acl_ {"ax_", "ay_", "ax_", "a_"};
  const std::vector<std::string> recognized_acl  {"ax",  "ay",  "ax",  "a"};
  const std::vector<std::string> recognized_vel  {"vx",  "vy",  "vx",  "v"};
  const std::vector<std::string> recognized_dsp  {"ux",  "uy",  "ux",  "u"};

  /* Advance all recognized fields. */
  for (PetscInt i = 0; i < recognized_acl.size(); i++) {
    if (fields.count(recognized_acl[i])) {
      VecAXPBYPCZ(fields[recognized_vel[i]]->mGlb, acl_factor, acl_factor, 1.0,
                  fields[recognized_acl[i]]->mGlb, fields[recognized_acl_[i]]->mGlb);
      VecAXPBYPCZ(fields[recognized_dsp[i]]->mGlb, mDt, dsp_factor, 1.0,
                  fields[recognized_vel[i]]->mGlb, fields[recognized_acl[i]]->mGlb);
      VecCopy(fields[recognized_acl[i]]->mGlb, fields[recognized_acl_[i]]->mGlb);
    } else { continue; }
  }

  time += mDt;
  return std::tuple<FieldDict, PetscScalar> (std::move(fields), time);

}
Order2Newmark::Order2Newmark(const std::unique_ptr<Options> &options) : Problem(options) {
  mDt = options->TimeStep();
}
