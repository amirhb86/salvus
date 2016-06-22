#include <Problem/ProblemNew.h>
#include <Utilities/Options.h>
#include <Utilities/Logging.h>
#include <Utilities/Utilities.h>
#include <Problem/Order2Newmark.h>
#include <Element/Element.h>
#include <Mesh/Mesh.h>
#include <stdexcept>

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

    LOG() << e.what();

  }

  return problem;

}

std::vector<std::unique_ptr<Element>> ProblemNew::initializeElements(std::unique_ptr<Mesh> &mesh,
                                                                     std::unique_ptr<ExodusModel> &model,
                                                                     std::unique_ptr<Options> const &options) {

  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  auto srcs = Source::Factory(options);
  auto recs = Receiver::Factory(options);
  std::vector<std::unique_ptr<Element>> elements;

  bool true_attach = true;
  bool trial_attach = false;
  std::vector<PetscInt> srcs_this_partition(srcs.size(), false);
  std::vector<PetscInt> recs_this_partition(recs.size(), false);

  /* Allocate all elements. */
  for (PetscInt i = 0; i < mesh->NumberElementsLocal(); i++)
  {

    /* Push back an appropriate element based on the mesh. */
    elements.push_back(Element::Factory({"u"}, {}, options));

    /* Assign a (processor-specific) number to this element. */
    elements.back()->SetNum(i);

    /* Attach vertex co-ordinates from mesh. */
    elements.back()->attachVertexCoordinates(mesh);

    /* Attach material properties (velocity, Cij, etc...). */
    elements.back()->attachMaterialProperties(model);

    /* Set any (external) boundary conditions. */
    elements.back()->setBoundaryConditions(mesh);

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
    std::cout << e.what() << std::endl;
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }

  return elements;

}

