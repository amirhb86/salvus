//
// Created by Michael Afanasiev on 2016-01-29.
//

#include "Utilities.h"

void utilities::print_from_root_mpi(const std::string msg) {
    if (MPI::COMM_WORLD.Get_rank() == 0) {
        std::cout << msg << std::endl;
    }
}

int utilities::broadcastInt(int send_buffer) {
    int root = 0;
    int num_ints = 1;
    MPI::COMM_WORLD.Bcast(&send_buffer, num_ints, MPI::INT, root);
    return send_buffer;
}

std::vector<double> utilities::broadcastStdVecFromRoot(std::vector<double> &send_buffer) {

    int root = 0;
    int int_size = 1;
    std::vector<double> receive_buffer;

    int length;
    if (!MPI::COMM_WORLD.Get_rank()) { length = send_buffer.size(); }
    MPI::COMM_WORLD.Bcast(&length, int_size, MPI::INT, root);

    receive_buffer.resize(length);
    if (!MPI::COMM_WORLD.Get_rank()) { receive_buffer = send_buffer; }
    MPI::COMM_WORLD.Bcast(receive_buffer.data(), length, MPI::DOUBLE, root);
    return receive_buffer;

}

std::vector<std::string> utilities::broadcastStringVecFromFroot(std::vector<std::string> &send_buffer) {

    int root = 0;
    int int_size = 1;

    // Find out how many string are stored in the vector.
    int length;
    if (!MPI::COMM_WORLD.Get_rank()) { length = send_buffer.size(); }
    MPI::COMM_WORLD.Bcast(&length, int_size, MPI::INT, root);

    // Broadcast each of these string individually.
    std::vector<std::string> receive_buffer;
    for (auto i = 0; i < length; i++) {

        // First broadcast the size of each string.
        int string_size;
        if (!MPI::COMM_WORLD.Get_rank()) { string_size = send_buffer[i].size(); }
        MPI::COMM_WORLD.Bcast(&string_size, int_size, MPI::INT, root);

        // Allocate receiving buffer to size + 1 for null c terminator.
        char *receive_c_buffer = new char [string_size + 1];
        if (!MPI::COMM_WORLD.Get_rank()) { send_buffer[i].copy(receive_c_buffer, string_size, 0); }
        MPI::COMM_WORLD.Bcast(receive_c_buffer, string_size, MPI::CHAR, root);

        // Add null c terminator to end of string and clean up.
        receive_c_buffer[string_size] = '\0';
        receive_buffer.push_back(std::string(receive_c_buffer));
        delete [] receive_c_buffer;
    }


    return receive_buffer;
}