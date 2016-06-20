#pragma once

// 3rd party.
#include <hdf5.h>

// parents.
#include <Receiver/Receiver.h>

// forward decl.
class Options;

class ReceiverHdf5 : public Receiver {

  static hid_t mFileId;
  hsize_t max_size;
  static std::vector<std::string> mWriteRegisteredFields;

public:

  ReceiverHdf5(std::unique_ptr<Options> const &options);
  ~ReceiverHdf5();

  void write();

};

