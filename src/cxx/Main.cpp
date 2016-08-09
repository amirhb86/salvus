#include <vector>
#include <memory>
#include <iostream>

#include <petsc.h>
#include <salvus.h>

INIT_LOGGING_STATE();

int main(int argc, char *argv[]) {

  try {

    /* Initialize PETSc, MPI, and command line args. */
    PetscInitialize(&argc, &argv, NULL, NULL);

    /* Parse command line options. */
    std::unique_ptr<Options> options(new Options);
    options->setOptions();

    /* Use options to allocate simulation components. */
    std::unique_ptr<Mesh> mesh(Mesh::Factory(options));
    std::unique_ptr<Problem> problem(Problem::Factory(options));
    std::unique_ptr<ExodusModel> model(new ExodusModel(options));

    /* Initialize relevant components and perform parallel decomposition. */
    mesh->read();
    model->read();

    /* Attach physics. Use this to inform the element generation. */
    mesh->setupTopology(model, options);

    /* Use mesh topology to generate our master list of elements. */
    auto elements = problem->initializeElements(mesh, model, options);

    /* Use elements to inform the global DOF layout. */
    mesh->setupGlobalDof(elements[0], options);

    auto fields = problem->initializeGlobalDofs(elements, mesh);

    /* Compute solution in time. */
    /* TODO: Make this something like problem->solve()? */
    PetscReal time = 0;
    while (time < options->Duration()) {

      /* Sum up all forces. */
      std::tie(elements, fields) = problem->assembleIntoGlobalDof(
          std::move(elements), std::move(fields), time,
          mesh->DistributedMesh(), mesh->MeshSection(), options);

      /* Apply inverse mass matrix. */
      fields = problem->applyInverseMassMatrix(std::move(fields));

      /* Advance time. */
      std::tie(fields, time) = problem->takeTimeStep(
          std::move(fields), time, options);

      problem->saveSolution(
          time, options->MovieFields(), fields, mesh->DistributedMesh());

      if (!PetscGlobalRank) { std::cout << "TIME: " << time << '\r'; std::cout.flush(); }

    }
  }

  /* TODO: Better MPI error handling. */
  catch (std::runtime_error &e) {
    LOG() << e.what();
    PetscFinalize();
    exit(1);
  }

  PetscFinalize();

}
