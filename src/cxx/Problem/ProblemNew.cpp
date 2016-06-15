#include <Problem/ProblemNew.h>
#include <Utilities/Options.h>
#include <Utilities/Logging.h>
#include <Problem/Order2Newmark.h>
#include <Element/Element.h>
#include <Mesh/Mesh.h>

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

  std::vector<std::unique_ptr<Element>> elements;

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


  }

  return elements;

}

