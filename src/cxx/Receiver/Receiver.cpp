#include "Receiver.h"
#include <iostream>

std::vector<Receiver *> Receiver::factory(Options options) {

  std::vector<Receiver *> receivers;
  for (int i = 0; i < options.NumberReceivers(); i++) {
    try {
      if (options.ReceiverType() == "hdf5") {

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

