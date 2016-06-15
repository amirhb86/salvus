#include <ctime>

#include <mpi.h>
#include <Mesh/Mesh.h>
#include <Source/Source.h>
#include <Element/Element.h>
#include <Receiver/Receiver.h>
#include <Model/ExodusModel.h>
#include <Utilities/Options.h>
#include <Utilities/Logging.h>
#include <Problem/NewmarkGeneral.h>

using namespace Eigen;

void NewmarkGeneral::initialize(Mesh *mesh,
                                std::unique_ptr<ExodusModel> const &model,
                                std::unique_ptr<Options> const &options) {

  // Save references to mesh and element base.
  mMesh = mesh;

  // Attach elements to mesh.
  if (options->Dimension() == 2) {
    mMesh->setupGlobalDof(1, 3, 9, 0, 2, model);
  } else if (options->Dimension() == 3) {
    mMesh->setupGlobalDof(1, 2, 4, 8, 3, model);
  }

  // Setup boundary conditions from options.
  mMesh->setupBoundaries(options);

  // Register all global field required for time stepping.
  for (auto field : mMesh->GlobalFields()) {
    mMesh->registerFieldVectors(field);
  }
  // Get a list of all local elements.

  for (PetscInt i = 0; i < mesh->NumberElementsLocal(); i++) {
    mElements.push_back(Element::Factory(mesh->ElementFields(i),
//                                         {}, options));
                                         mesh->TotalCouplingFields(i),
                                         options));
  }


  // Get a list of all sources and receivers.
  auto sources = Source::factory(options);
  mRecs = Receiver::factory(options);

  // Set up elements.
  Eigen::VectorXd h_all(mElements.size());
  int element_number = 0;
  for (auto &element : mElements) {

    // Give each element a number starting from zero.
    element->SetNum(element_number); element_number++;

    // Get vertex coordinates from the PETSc DMPLEX.
    element->attachVertexCoordinates(mMesh);

    // Set boundary conditions.
    element->setBoundaryConditions(mMesh);

    // Add material parameters (velocity, Cij, etc...).
    element->attachMaterialProperties(model);

    // Assemble the (elemental) mass matrix.
    element->assembleElementMassMatrix(mMesh);

    // Attach any external sources and receivers.
    element->attachSource(sources);
    element->attachReceiver(mRecs);

    // Prepare stiffness terms.
    element->prepareStiffness();

    // Save CFL estimate.
    h_all(element_number - 1) = element->CFL_estimate();

  }
  
  /**** Time step ****/
  auto dt = mesh->CFL() * h_all.minCoeff();
  // if timestep not set from command line
  if(options->TimeStep() <= 0) {
//    options.SetTimeStep(dt);
    LOG() << "Suggested dt = " << options->TimeStep();
  } else {
    LOG() << "Suggested dt = " << dt << " vs. Commandline dt = " << options->TimeStep();
  }
  /**** END NEW TIME STEP ****/

  // Scatter the mass matrix to the global dofs.
  mMesh->checkInFieldBegin("m");
  mMesh->checkInFieldEnd("m");

  // TODO: NOTE HERE! Need to make this a smart pointer as well.

}

void NewmarkGeneral::solve(Options options) {

  std::set<std::string> push_fields, pull_fields;
  if (options.MeshType() == "3d_couple") {
    push_fields = {"a", "ax", "ay", "az"};
    pull_fields = {"vx", "vy", "vz", "v", "ux", "uy", "uz", "u"};
  } else if (options.MeshType() == "2d_couple") {
    push_fields = {"a", "ax", "ay"};
    pull_fields = {"vx", "vy", "v", "ux", "uy", "u"};
  }

  // Setup values.
  int it = 0;
  int it_movie = 0;
  double time = 0.0;
  double timeStep = options.TimeStep();
  double duration = options.Duration();
  if (options.SaveMovie()) mMesh->setUpMovie(options.OutputMovieFile());

  // Over-allocate matrices to avoid re-allocations.
  int max_dims = 3;
  int int_pnts =  mElements.front()->NumIntPnt();
  Eigen::MatrixXd s(int_pnts, max_dims);
  Eigen::MatrixXd f(int_pnts, max_dims);
  Eigen::MatrixXd u(int_pnts, max_dims);
  Eigen::MatrixXd ku(int_pnts, max_dims);
  Eigen::MatrixXd fMinusKu(int_pnts, max_dims);
  double max_Loo = 0.0;

  clock_t time_start = clock();
  while (time < duration) {
    max_Loo = 0.0;

    // Collect all global fields to the local partitions.
    for (auto &field : pull_fields) {
      mMesh->checkOutField(field);
    }

    // Zero all fields to which we will sum.
    for (auto &field : push_fields) {
      mMesh->zeroFields(field);
    }

    for (auto &element : mElements) {

      /* The number of fields to pull may be different from those you push.
       * This is the case on coupling interfaces. */
      int num_pull_fields = element->PullElementalFields().size();
      int num_push_fields = element->PushElementalFields().size();

      // Get fields on element, store in successive rows.
      int fitr = 0;
      for (auto &field : element->PullElementalFields()) {
        u.col(fitr) = mMesh->getFieldOnElement(field, element->Num(),
                                               element->ClsMap());
        fitr++;
        if (options.DisplayDiagnostics() && (it % options.DisplayDiagnosticsEvery() == 0 || it == 0)) {
          max_Loo = fmax(max_Loo, u.col(fitr - 1).array().abs().maxCoeff());
        }
      }

      // Record displacement.
      element->recordField(u.block(0, 0, int_pnts, num_pull_fields));

      // Compute stiffness, only passing those rows which are occupied.
      ku.leftCols(num_push_fields) =
          element->computeStiffnessTerm(u.leftCols(num_pull_fields));


      // Compute source term.
      f.leftCols(num_push_fields) = element->computeSourceTerm(time);

      // Compute surface integral.
      s.leftCols(num_push_fields) = element->computeSurfaceIntegral(u.leftCols(num_pull_fields));

      // Compute acceleration.
      fMinusKu.leftCols(num_push_fields) = f.leftCols(num_push_fields).array() -
          ku.leftCols(num_push_fields).array() + s.leftCols(num_push_fields).array();

      // Sum fields into local partition.
      fitr = 0;
      for (auto &field : element->PushElementalFields()) {
        mMesh->addFieldFromElement(field, element->Num(),
                                   element->ClsMap(),
                                   fMinusKu.col(fitr));
        fitr++;
      }
    }

//    if (it > 0)
//    exit(0);


    PetscInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (options.DisplayDiagnostics() && (it % options.DisplayDiagnosticsEvery() == 0 || it == 0)) {
      //  if(rank == 0) printf("|u|_oo=%f @ time=%f (%f%%)\n",max_Loo,time,100*(time/duration));
    }

    // Sum fields into global partitions.
    for (auto &field : push_fields) {
      mMesh->checkInFieldBegin(field);
      mMesh->checkInFieldEnd(field);
    }

    // Take a time step.
    mMesh->applyInverseMassMatrix();
    mMesh->advanceField(timeStep);

    if (options.SaveMovie() && (it % options.SaveFrameEvery() == 0 || it == 0)) {
      // GlobalFields[0] == "u" for acoustic and == "ux" for elastic
      mMesh->saveFrame("u", it_movie);
      mMesh->saveFrame("ux", it_movie);
      it_movie++;
    }

    it++;
    time += timeStep;
    if (rank == 0) std::cout << "TIME: " << time << "\n" << std::flush;

  }

  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank)
    printf("Time elapsed: %f\n", ((double)clock() - time_start) / CLOCKS_PER_SEC);

  if (options.SaveMovie()) mMesh->finalizeMovie();

  // Save all receivers.
  for (auto &rec : mRecs) {
    rec->write();
  }

}

