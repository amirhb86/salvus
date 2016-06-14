#include <Problem/ProblemNew.h>
#include <Utilities/Options.h>
#include <Utilities/Logging.h>
#include <Problem/Order2Newmark.h>
#include <Element/Element.h>

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
                                                                     std::unique_ptr<Model> &model,
                                                                     std::unique_ptr<Options> const &options) {

  std::vector<std::unique_ptr<Element>> elements;
  for (PetscInt i = 0; i < 10; i++) {
    elements.push_back(Element::Factory({"u"}, {}, options));
    elements.back()->SetNum(i);
  }
  return elements;

}

