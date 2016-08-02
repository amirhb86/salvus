#include <iostream>
#include <assert.h>
#include <Utilities/Options.h>
#include <Receiver/Receiver.h>
#include <Receiver/ReceiverHdf5.h>
#include <stdexcept>
#include <Utilities/Utilities.h>
#include <Utilities/Logging.h>

// Initialize counter to zero.
PetscInt Receiver::mNumRecs = 0;

std::vector<std::unique_ptr<Receiver>> Receiver::Factory(std::unique_ptr<Options> const &options) {

  std::vector<std::unique_ptr<Receiver>> receivers;
  for (int i = 0; i < options->NumberReceivers(); i++) {
    if (utilities::stringHasExtension(options->ReceiverFileName(), ".h5")) {
      receivers.push_back(std::unique_ptr<ReceiverHdf5>(new ReceiverHdf5(options)));
    } else {
      throw std::runtime_error("Runtime error: Filetype of receiver file cannot be deduced from extension."
                                   " Use [ .h5 ]");
    }
  }
  return receivers;
}

Receiver::Receiver(std::unique_ptr<Options> const &options) {

  // Get receiver number and increment.
  SetNum(mNumRecs++);

  // Set physical location.
  mLocX = options->RecLocX()[mNum];
  mLocY = options->RecLocY()[mNum];
  if (options->RecLocZ().size()) {
    mLocZ = options->RecLocZ()[mNum];
  }

  // Set name.
  mName = options->RecNames()[mNum];

}



Receiver::~Receiver() { --mNumRecs; }

void Receiver::record(const double val, const std::string &field) {

  // If entry does not exist, create it.
  if (store.find(field) == store.end()) {
    store.insert(std::pair<std::string, std::vector<float>> (field, std::vector<float> {}));
  }

  // Push back to vector.
  store[field].push_back(val);

}

