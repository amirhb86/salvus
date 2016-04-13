#include <petscsys.h>
#include "Utilities.h"

void utilities::print_from_root_mpi(const std::string msg) {
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << msg << std::endl;
  }
}

utilities::PrintFromRoot::PrintFromRoot() { }

utilities::PrintFromRoot::~PrintFromRoot() {
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << os.str() << std::endl;
  }
}

int utilities::broadcastInt(int send_buffer) {
  int root = 0;
  int num_ints = 1;
  MPI_Bcast(&send_buffer, num_ints, MPI_INT, root, PETSC_COMM_WORLD);
  return send_buffer;
}

std::vector<double> utilities::broadcastStdVecFromRoot(std::vector<double> &send_buffer) {

  int root = 0;
  int int_size = 1;
  std::vector<double> receive_buffer;

  long length;
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) { length = send_buffer.size(); }
  MPI_Bcast(&length, int_size, MPI_LONG, root, PETSC_COMM_WORLD);

  receive_buffer.resize(length);
  if (!rank) { receive_buffer = send_buffer; }
  MPI_Bcast(receive_buffer.data(), length, MPI_DOUBLE, root, PETSC_COMM_WORLD);
  return receive_buffer;

}

std::vector<int> utilities::broadcastStdVecFromRoot(std::vector<int> &send_buffer) {

  int root = 0;
  int int_size = 1;
  std::vector<int> receive_buffer;

  long length;
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) { length = send_buffer.size(); }
  MPI_Bcast(&length, int_size, MPI_LONG, root, PETSC_COMM_WORLD);

  receive_buffer.resize(length);
  if (!rank) { receive_buffer = send_buffer; }
  MPI_Bcast(receive_buffer.data(), length, MPI_INT, root, PETSC_COMM_WORLD);
  return receive_buffer;

}

std::vector<std::string> utilities::broadcastStringVecFromRank(std::vector<std::string> &send_buffer,
                                                               int root) {

  int int_size = 1;
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Find out how many string are stored in the vector.
  long length;
  if (rank == root) { length = send_buffer.size(); }
  MPI_Bcast(&length, int_size, MPI_LONG, root, PETSC_COMM_WORLD);

  // Broadcast each of these string individually.
  std::vector<std::string> receive_buffer;
  for (auto i = 0; i < length; i++) {

    // First broadcast the size of each string.
    long string_size;
    if (rank == root) { string_size = send_buffer[i].size(); }
    MPI_Bcast(&string_size, int_size, MPI_LONG, root, PETSC_COMM_WORLD);

    // Allocate receiving buffer to size + 1 for null c terminator.
    char *receive_c_buffer = new char[string_size + 1];
    if (rank == root) { send_buffer[i].copy(receive_c_buffer, string_size, 0); }
    MPI_Bcast(receive_c_buffer, string_size, MPI_CHAR, root, PETSC_COMM_WORLD);

    // Add null c terminator to end of string and clean up.
    receive_c_buffer[string_size] = '\0';
    receive_buffer.push_back(std::string(receive_c_buffer));
    delete[] receive_c_buffer;
  }

  return receive_buffer;
}

std::string utilities::broadcastStringFromRank(std::string &buf, int rank) {

  // Get size of string from rank.
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  long name_size = (rank == my_rank) ? buf.size() : 0;
  MPI_Bcast(&name_size, 1, MPI_LONG, rank, PETSC_COMM_WORLD);

  char *receive_c_buffer = new char[name_size+1];
  if (my_rank == rank) { buf.copy(receive_c_buffer, name_size, 0); }
  MPI_Bcast(receive_c_buffer, name_size, MPI_CHAR, rank, PETSC_COMM_WORLD);

  // Add null terminator.
  receive_c_buffer[name_size] = '\0';
  std::string ret_string = std::string(receive_c_buffer);

  delete [] receive_c_buffer;
  return ret_string;


}


