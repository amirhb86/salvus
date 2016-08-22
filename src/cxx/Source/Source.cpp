#include <mpi.h>
#include <Source/Source.h>
#include <Source/Ricker.h>
#include <Source/SourceHdf5.h>
#include <Utilities/Options.h>
#include <stdexcept>
#include <iostream>
#include <Utilities/Logging.h>

/* Initialize counter. */
PetscInt Source::number = 0;

enum source_type { sRicker, sHDF5, sTypeError };
source_type stype(const std::string &stype) {
  if (stype == "ricker") return sRicker;
  if (stype == "file")  return sHDF5;
  return sTypeError;
}

std::vector<std::unique_ptr<Source>> Source::Factory(std::unique_ptr<Options> const &options) {

  std::vector<std::unique_ptr<Source>> sources;
  if ( options->NumberSources() > 0) {
    switch (stype(options->SourceType())) {

      case sRicker:
        for (PetscInt i = 0; i < options->NumberSources(); i++) {
          sources.push_back(std::unique_ptr<Ricker>(new Ricker(options)));
        }
        return sources;

      case sHDF5:
        for (PetscInt i = 0; i < options->NumberSources(); i++) {
          sources.push_back(std::unique_ptr<SourceHdf5>(new SourceHdf5(options)));
        }
        return sources;

      case sTypeError:
        break;
    }

    throw std::runtime_error("Source type '" + options->SourceType() + "' not supported.");
  }
  return sources;
}

Source::Source(std::unique_ptr<Options> const &options) {

  /* Increment global number, save this particular one. */
  SetNum(number++);

}

Source::~Source() { --number; }

void Source::loadData() {}
