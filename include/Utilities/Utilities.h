#pragma once

// stl.
#include <iosfwd>
#include <string>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

// 3rd party.
#include <mpi.h>

namespace utilities {

void print_from_root_mpi(const std::string msg);

int broadcastInt(int send_buffer);
std::vector<int> broadcastStdVecFromRoot(std::vector<int> &send_buffer);
std::vector<double> broadcastStdVecFromRoot(std::vector<double> &send_buffer);
std::string broadcastStringFromRank(std::string &buf, int rank);
std::vector<std::string> broadcastStringVecFromRank(std::vector<std::string> &send_buffer, int rank);


class PrintFromRoot {

 protected:
  std::ostringstream os;
 private:
  PrintFromRoot(const PrintFromRoot &);
  PrintFromRoot &operator=(const PrintFromRoot &);

 public:
  PrintFromRoot();
  virtual ~PrintFromRoot();
  std::ostringstream &Get() { return os; }
};

}

#define PRINT_ROOT() utilities::PrintFromRoot().Get()

enum class AcousticFields {
  displacement = 0, velocity = 1, acceleration = 2, force = 3, acceleration_ = 4, mass_matrix = -1,
  mass_matrix_inverse = -2
};
