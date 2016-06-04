#include <mpi.h>
#include <Mesh/Mesh.h>
#include <Source/Source.h>
#include <Element/Element.h>
#include <Receiver/Receiver.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Problem/NewmarkGeneral.h>

using namespace Eigen;

void NewmarkGeneral::initialize(Mesh *mesh,
                                ExodusModel *model,
                                std::shared_ptr<Element> elem,
                                Options options) {

  // Save references to mesh and element base.
  mMesh = mesh;
  mReferenceElem = elem;
  // Attach elements to mesh.
  mMesh->setupGlobalDof(mReferenceElem->NumDofVtx(),
                        mReferenceElem->NumDofEdg(),
                        mReferenceElem->NumDofFac(),
                        mReferenceElem->NumDofVol(),
                        mReferenceElem->NumDim(),
                        model);

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

  // Get a list of all sources and receivers.
  auto sources = Source::factory(options);
  mRecs = Receiver::factory(options);

  // Set up elements.
  int element_number = 0;
  for (auto &element : mElements) {

    // Give each element a number starting from zero.
    element->SetNum(element_number); element_number++;

    // Get vertex coordinates from the PETSc DMPLEX.
    element->attachVertexCoordinates(mMesh);

    // Add material parameters (velocity, Cij, etc...).
    element->attachMaterialProperties(model);

    // Set boundary conditions.
    element->setBoundaryConditions(mMesh);

    // Assemble the (elemental) mass matrix.
    element->assembleElementMassMatrix(mMesh);

    // Attach any external sources and receivers.
    element->attachSource(sources);
    element->attachReceiver(mRecs);

    // Prepare stiffness terms.
    element->prepareStiffness();

  }

  // Scatter the mass matrix to the global dofs.
  mMesh->checkInFieldBegin("m");
  mMesh->checkInFieldEnd("m");

  // TODO: NOTE HERE! Need to make this a smart pointer as well.
  delete model;

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
  Eigen::MatrixXd s(int_pnts, max_dims);
  Eigen::MatrixXd f(int_pnts, max_dims);
  Eigen::MatrixXd u(int_pnts, max_dims);
  Eigen::MatrixXd ku(int_pnts, max_dims);
  Eigen::MatrixXd fMinusKu(int_pnts, max_dims);
  double max_Loo = 0.0;
  while (time < duration) {
    max_Loo = 0.0;
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
        if(options.DisplayDiagnostics() && (it % options.DisplayDiagnosticsEvery() == 0 || it == 0)) {
          max_Loo = fmax(max_Loo,u.col(fitr-1).array().abs().maxCoeff());
        }
      }

      // Record displacement.
      element->recordField(u.block(0, 0, int_pnts, fitr));

      // Compute stiffness, only passing those rows which are occupied.
      ku.leftCols(fitr) =
        element->computeStiffnessTerm(u.leftCols(fitr));

      // Compute source term.
      f.leftCols(fitr) = element->computeSourceTerm(time);

      // Compute surface integral.
      s.leftCols(fitr) = element->computeSurfaceIntegral(u.leftCols(fitr));

      // Compute acceleration.
      fMinusKu.leftCols(fitr) = f.leftCols(fitr).array() -
        ku.leftCols(fitr).array();

      // Sum fields into local partition.
      fitr = 0;
      for (auto &field : mReferenceElem->PushElementalFields()) {
        mMesh->addFieldFromElement(field, element->Num(),
                                   element->ClsMap(),
                                   fMinusKu.col(fitr));
        fitr++;
      }
    }

    PetscInt rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if(options.DisplayDiagnostics() && (it % options.DisplayDiagnosticsEvery() == 0 || it == 0)) {
      if(rank == 0) printf("|u|_oo=%f @ time=%f (%f%%)\n",max_Loo,time,100*(time/duration));
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
      if(MPI::COMM_WORLD.Get_rank() == 0) printf("Saving frame %d\n",it);
      mMesh->saveFrame(mMesh->GlobalFields()[0], it);
    }

    it++;
    time += timeStep;
  }

  if (options.SaveMovie()) mMesh->finalizeMovie();

  // Save all receivers.
  for (auto &rec : mRecs) {
    rec->write();
  }

}

