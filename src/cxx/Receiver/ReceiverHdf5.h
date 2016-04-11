#pragma once

#include <hdf5.h>
#include <Utilities/Options.h>
#include "Receiver.h"
#include "ReceiverHdf5.h"

class ReceiverHdf5 : public Receiver {

  static hid_t mPlistId, mFileId;

 public:
  ReceiverHdf5(Options options);

};

