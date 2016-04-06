//
// Created by Michael Afanasiev on 2016-03-14.
//

#include "NewmarkGeneral.h"

using namespace Eigen;

void NewmarkGeneral::initialize(Mesh *mesh,
                                ExodusModel *model,
                                Element *elem,
                                Options options) {

  // Save references to mesh and element base.
  mMesh = mesh;
  mReferenceElem = elem;

  // Attach elements to mesh.
  mMesh->setupGlobalDof(mReferenceElem->NumDofVtx(),
                        mReferenceElem->NumDofEdg(),
                        mReferenceElem->NumDofFac(),
                        mReferenceElem->NumDofVol(),
                        mReferenceElem->NumDim());

  // Setup boundary conditions from options.
  mMesh->setupBoundaries(options);

  // Register all global field required for time stepping.
  for (auto field : mMesh->GlobalFields()) {
    mMesh->registerFieldVectors(field);
  }

  // Get a list of all local elements.
  for (int i = 0; i < mMesh->NumberElementsLocal(); i++) {
    mElements.push_back(mReferenceElem->clone());
  }

  // Get a list of all sources.
  auto sources = Source::factory(options);

  // Set up elements.
  int element_number = 0;
  for (auto &element : mElements) {

    // Give each element a number starting from zero.
    element->SetNum(element_number++);

    // Get vertex coordinates from the PETSc DMPLEX.
    element->attachVertexCoordinates(mMesh->DistributedMesh());

    // Add material parameters (velocity, Cij, etc...).
    element->interpolateMaterialProperties(model);

    // Set boundary conditions.
    element->setBoundaryConditions(mMesh);

    // Assemble the (elemental) mass matrix.
    element->assembleElementMassMatrix(mMesh);

    // Attach any external source terms.
    element->attachSource(sources);

    // Prepare stiffness terms.
    element->prepareStiffness();

  }

  // Scatter the mass matrix to the global dofs.
  mMesh->checkInFieldBegin("m");
  mMesh->checkInFieldEnd("m");

}

void NewmarkGeneral::solve(Options options) {

  // Setup values.
  int it = 0;
  double time = 0.0;
  double timeStep = options.TimeStep();
  double duration = options.Duration();
  if (options.SaveMovie()) mMesh->setUpMovie(options.OutputMovieFile());

  // Over-allocate matrices to avoid re-allocations.
  int max_dims = 3;
  int int_pnts = mReferenceElem->NumIntPnt();
  Eigen::MatrixXd f(int_pnts, max_dims);
  Eigen::MatrixXd u(int_pnts, max_dims);
  Eigen::MatrixXd ku(int_pnts, max_dims);
  Eigen::MatrixXd fMinusKu(int_pnts, max_dims);

  while (time < duration) {

    // Collect all global fields to the local partitions.
    for (auto &field : mReferenceElem->PullElementalFields()) {
      mMesh->checkOutField(field);
    }

    // Zero all fields to which we will sum.
    for (auto &field : mReferenceElem->PushElementalFields()) {
      mMesh->zeroFields(field);
    }

    for (auto &element : mElements) {

      // Get fields on element, store in successive rows.
      int fitr = 0;
      for (auto &field : mReferenceElem->PullElementalFields()) {
        u.col(fitr) = mMesh->getFieldOnElement(field, element->Num(),
                                               element->ClsMap());
        fitr++;
      }

      // Compute stiffness, only passing those rows which are occupied.
      ku.block(0, 0, int_pnts, fitr) = element->computeStiffnessTerm(
          u.block(0, 0, int_pnts, fitr));

      // Compute source term.
      f.block(0, 0, int_pnts, fitr) = element->computeSourceTerm(time);

      // Compute acceleration.
      fMinusKu.block(0, 0, int_pnts, fitr) = f.block(0, 0, int_pnts, fitr).array() -
          ku.block(0, 0, int_pnts, fitr).array();

      // Sum fields into local partition.
      fitr = 0;
      for (auto &field : mReferenceElem->PushElementalFields()) {
        mMesh->addFieldFromElement(field, element->Num(),
                                   element->ClsMap(),
                                   fMinusKu.col(fitr));
        fitr++;
      }
    }

    // Sum fields into global partitions.
    for (auto &field : mReferenceElem->PushElementalFields()) {
      mMesh->checkInFieldBegin(field);
      mMesh->checkInFieldEnd(field);
    }

    // Take a time step.
    mMesh->applyInverseMassMatrix();
    mMesh->advanceField(timeStep);

    if (options.SaveMovie() && (it % options.SaveFrameEvery() == 0 || it == 0)) {
      // GlobalFields[0] == "u" for acoustic and == "ux" for elastic
      mMesh->saveFrame(mMesh->GlobalFields()[0], it);
    }

    it++;
    time += timeStep;
    PRINT_ROOT() << "TIME: " << time;
  }

  if (options.SaveMovie()) mMesh->finalizeMovie();

}
