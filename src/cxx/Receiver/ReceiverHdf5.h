#pragma once

#include <hdf5.h>
#include <Utilities/Options.h>
#include "Receiver.h"
#include "ReceiverHdf5.h"

class ReceiverHdf5 : public Receiver {

  static hid_t mFileId;
  hsize_t max_size;
  static std::vector<std::string> mWriteRegisteredFields;

public:

  std::shared_ptr<Receiver> clone() const {
    return (std::shared_ptr<Receiver> (new ReceiverHdf5(*this)));
  }
  ReceiverHdf5(Options options);
  ~ReceiverHdf5();

  void write();

};

