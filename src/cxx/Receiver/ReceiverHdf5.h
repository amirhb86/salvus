#pragma once

#include <Utilities/Options.h>
#include "Receiver.h"
#include "ReceiverHdf5.h"

class ReceiverHdf5 : public Receiver {

 public:

  ReceiverHdf5(Options options);

};

