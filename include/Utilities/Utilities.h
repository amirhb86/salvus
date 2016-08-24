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
#include <Eigen/Dense>
#include <set>
#include <Element/Element.h>
#include "Types.h"

namespace utilities {

  /**
   * Database of (tag, {ranks}) pairs. Used to identify which ranks a certain mixin resides on.
   */
  extern std::set<std::tuple<std::string, std::vector<PetscInt>>> global_parallel_tags;

  /**
   * Small function to return whenter or not a string has a certain suffix.
   * @param str String to check (i.e. filename).
   * @param ext Wanted extension.
   * @return True if extension exists, false if not.
   */
  bool stringHasExtension(const std::string &str, const std::string &ext);

  /**
   * Small function to return an MPI datatype based on a template datatype.
   */
  template <typename T>
  MPI_Datatype mpitype();

  /**
   * Broadcasts a number on some rank to all processors.
   * @param [in] send_buffer Number to send.
   * @param [in] rank Rank to send from.
   * @returns Number broadcasted to all processors.
   */
  template <typename T>
  T broadcastNumberFromRank(T send_buffer, int rank);

  /**
  * Broadcasts a vector of numbers on some rank to all processors.
  * @param [in] send_buffer Vector of numbers to send.
  * @param [in] rank Rank to send from.
  * @returns Vector broadcasted to all processors.
  */
  template <typename T>
  std::vector<T> broadcastNumberVecFromRank(std::vector<T> &send_buffer, int rank);

  /**
   * c++ strings are a bit harder to broadcast than fundamental datatypes. This function
   * serializes a string, and broadcasts it to all ranks.
   * @param [in] buf String to send.
   * @param [in] rank Rank to send from.
   * @returns String broadcasted to all processors.
   */
  std::string broadcastStringFromRank(std::string &buf, int rank);

  /**
   * c++ string are a bit harder to broadcast than fundamental datatypes. This function
   * serializes a vector of strings, and broadcasts it to all ranks.
   * @param [in] send_buffer Vector of strings to send.
   * @param [in] rank Rank to send from.
   * @returns Vector of strings broadcasted to all processors.
   */
  std::vector<std::string> broadcastStringVecFromRank(std::vector<std::string> &send_buffer, int rank);

  /**
   * Given some tag (i.e. "SaveSurface"), returns the MPI ranks which contain elements containing
   * the tag in their name (SaveSurface_TensorQuad_QuadP1).
   * @param elements List of all elements.
   * @param tag Identifying tag.
   * @return Integer vector of ranks.
   */
  std::vector<PetscInt> getRanksForMixin(ElemVec &elements, const std::string &tag);

  /**
   * Insert an integer vector of MPI ranks into a global database which can be queried with
   * GetGlobalRankTags();
   * @param ranks Vector of ranks on which elements with `tag` reside.
   * @param tag Identifying tag.
   */
  void appendToGlobalRankTags(const std::vector<PetscInt> ranks, const std::string &tag);

  /**
   * Query the global database for the MPI ranks corresponding to a particular tag.
   * @param tag Identifyihng tag.
   * @return Integer vector of ranks corresponding to tag.
   */
  std::vector<PetscInt> GetWorldRanksForTag(const std::string &tag);


}

void seg_scan(int *invec, int *inoutvec, int *len, MPI_Datatype *dtype);

