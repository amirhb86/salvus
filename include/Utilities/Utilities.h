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

namespace utilities {

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

}

void seg_scan(int *invec, int *inoutvec, int *len, MPI_Datatype *dtype);

