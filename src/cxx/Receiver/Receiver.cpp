#include <iostream>
#include <assert.h>
#include <Utilities/Options.h>
#include <Receiver/Receiver.h>
#include <Receiver/ReceiverHdf5.h>
#include <stdexcept>

// Initialize counter to zero.
long Receiver::num = 0;

std::vector<std::shared_ptr<Receiver>> Receiver::factory(std::unique_ptr<Options> const &options) {

  std::vector<std::shared_ptr<Receiver>> receivers;
  for (int i = 0; i < options->NumberReceivers(); i++) {
    try {
      if (options->ReceiverType() == "hdf5") {
        receivers.push_back(std::shared_ptr<ReceiverHdf5>(new ReceiverHdf5(options)));
      } else {
        throw std::runtime_error("Runtime error: Receiver type " + options->ReceiverType() + " not supported.");
      }

    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  return receivers;
}

Receiver::Receiver(std::unique_ptr<Options> const &options) {

  // Check sanity of receiver coordinates.
  if (options->Dimension() == 2) {
    assert(options->RecLocX3().size() == 0);
  } else if (options->Dimension() == 3) {
    assert(options->RecLocX3().size() == options->RecLocX2().size());
  }


  // Get receiver number and increment.
  mNum = Receiver::num;
  Receiver::num++;

  // Set physical location.
  mPysLocX1 = options->RecLocX1()[mNum];
  mPysLocX2 = options->RecLocX2()[mNum];
  if (options->RecLocX3().size()) {
    mPysLocX3 = options->RecLocX3()[mNum];
  }

  // Set name.
  mName = options->RecNames()[mNum];

}



Receiver::~Receiver() {

  // Reset receiver count.
  Receiver::num = 0;

}
void Receiver::record(const double val, const std::string &field) {

  // If entry does not exist, create it.
  if (store.find(field) == store.end()) {
    store.insert(std::pair<std::string, std::vector<float>> (field, std::vector<float> {}));
  }

  // Push back to vector.
  store[field].push_back(val);

}

