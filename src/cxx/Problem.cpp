//
// Created by Michael Afanasiev on 2016-01-27.
//

#include <ostream>
#include <iostream>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include "Problem.h"
#include "Utilities.h"


Problem *Problem::factory(std::string solver_type) {
    try {

        if (solver_type == "time_domain") {
            return new TimeDomain;
        } else if (solver_type == "frequency_domain") {
            return new FrequencyDomain;
        } else {
            throw std::runtime_error("Runtime Error: Problem type " + solver_type + " not supported.");
        }
    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI::COMM_WORLD.Abort(-1);
        return nullptr;
    }
}

void FrequencyDomain::initialize() {

    utilities::print_from_root_mpi("Initializing Frequency Domain Problem.");
}

void FrequencyDomain::solve() { }

void TimeDomain::initialize() {

    utilities::print_from_root_mpi("Initializing Frequency Domain Problem.");

}

void TimeDomain::solve() { }