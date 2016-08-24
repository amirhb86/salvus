#include <mpi.h>
#include <petsc.h>
#include <Utilities/Utilities.h>
#include <Utilities/Types.h>

/* Template specifications for different MPI datatypes. */
namespace utilities {
  template<>
  MPI_Datatype mpitype<double>() { return MPI_DOUBLE; }

  template<>
  MPI_Datatype mpitype<float>() { return MPI_FLOAT; }

  template<>
  MPI_Datatype mpitype<PetscInt>() { return MPIU_INT; }

  std::set<std::tuple<std::string, std::vector<PetscInt>>> global_parallel_tags;

}

template <typename T>
T utilities::broadcastNumberFromRank(T send_buffer, int rank) {
  int num_ints = 1;
  MPI_Bcast(&send_buffer, num_ints, mpitype<T>(), rank, PETSC_COMM_WORLD);
  return send_buffer;
}

template <typename T>
std::vector<T> utilities::broadcastNumberVecFromRank(std::vector<T> &send_buffer, int rank) {

  int int_size = 1;
  std::vector<T> receive_buffer;

  long length;
  int myrank; MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
  if (myrank == rank) { length = send_buffer.size(); }
  MPI_Bcast(&length, int_size, MPI_LONG, rank, PETSC_COMM_WORLD);

  receive_buffer.resize(length);
  if (myrank == rank) { receive_buffer = send_buffer; }
  MPI_Bcast(receive_buffer.data(), length, mpitype<T>(), rank, PETSC_COMM_WORLD);
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

bool ::utilities::stringHasExtension(const std::string &str, const std::string &ext) {
  return str.size() >= ext.size() &&
      str.compare(str.size() - ext.size(), ext.size(), ext) == 0;
}

std::vector<PetscInt> utilities::getRanksForMixin(ElemVec &elements,
                                                  const std::string &tag) {

  PetscInt size; MPI_Comm_size(PETSC_COMM_WORLD, &size);
  PetscInt rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  std::vector<PetscInt> ranks(size, 0);
  for (auto &elm: elements) {
    std::string name = elm->Name();
    if (name.find(tag) != std::string::npos) { ranks[rank] = 1; }
  }
  MPI_Allreduce(MPI_IN_PLACE, ranks.data(), ranks.size(), MPIU_INT, MPI_MAX,
                PETSC_COMM_WORLD);
  PetscInt num = std::count_if(ranks.begin(), ranks.end(), [](PetscInt i) { return i == 1; });
  PetscInt j = 0, k = 0; std::vector<PetscInt> found_ranks(num);
  std::for_each(ranks.begin(), ranks.end(),
                [&found_ranks, &j, &k]
                    (PetscInt i) { if (i) found_ranks[k++] = j; j++; });
  return found_ranks;

}

void utilities::appendToGlobalRankTags(const std::vector<PetscInt> ranks, const std::string &tag) {
  utilities::global_parallel_tags.insert(
      std::tuple<std::string, std::vector<PetscInt>>(tag, ranks));
}

std::vector<PetscInt> utilities::GetWorldRanksForTag(const std::string &tag) {
  for (auto &tags: global_parallel_tags) {
    if (std::get<0>(tags) == tag) { return std::get<1>(tags); }
  }
  throw std::runtime_error("Tag " + tag + " not found in global database.");
}

template PetscInt utilities::broadcastNumberFromRank(PetscInt send_buffer, int rank);
template float utilities::broadcastNumberFromRank(float send_buffer, int rank);
template double utilities::broadcastNumberFromRank(double send_buffer, int rank);

template std::vector<PetscInt> utilities::broadcastNumberVecFromRank(std::vector<PetscInt> &send_buffer, int rank);
template std::vector<float> utilities::broadcastNumberVecFromRank(std::vector<float> &send_buffer, int rank);
template std::vector<double> utilities::broadcastNumberVecFromRank(std::vector<double> &send_buffer, int rank);


