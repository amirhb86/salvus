#include "Receiver.h"
#include <iostream>
#include <exception>

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
Receiver::Receiver(Options options) {

  mPysLocX = 0.0;
  mPysLocY = 0.0;
  mPysLocZ = 0.0;

}



