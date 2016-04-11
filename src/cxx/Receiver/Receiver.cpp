#include "Receiver.h"
#include "ReceiverHdf5.h"
#include <iostream>
#include <assert.h>

// Initialize counter to zero.
long Receiver::num = 0;

std::vector<std::shared_ptr<Receiver>> Receiver::factory(Options options) {

  std::vector<std::shared_ptr<Receiver>> receivers;
  for (int i = 0; i < options.NumberReceivers(); i++) {
    try {
      if (options.ReceiverType() == "hdf5") {
        receivers.push_back(std::shared_ptr<ReceiverHdf5> (new ReceiverHdf5(options)));
      } else {
        throw std::runtime_error("Runtime error: Receiver type " + options.ReceiverType() + " not supported.");
      }

    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

  }
  return receivers;
}

Receiver::Receiver(Options options) {

  // Check sanity of receiver coordinates.
  if (options.Dimension() == 2) {
    assert(options.RecLocX3().size() == 0);
  } else if (options.Dimension() == 3) {
    assert(options.RecLocX3().size() == options.RecLocX2().size());
  }

  // Get receiver number and increment.
  mNum = Receiver::num;
  Receiver::num++;

  // Set physical location.
  mPysLocX1 = options.RecLocX1()[mNum];
  mPysLocX2 = options.RecLocX2()[mNum];
  if (options.RecLocX3().size()) {
    mPysLocX3 = options.RecLocX3()[mNum];
  }

}
Receiver::~Receiver() {

  // Reset receiver count.
  Receiver::num = 0;

}





