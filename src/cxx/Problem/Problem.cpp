//
// Created by Michael Afanasiev on 2016-01-27.
//

#include "Problem.h"
#include "NewmarkGeneral.h"
#include "NewmarkTesting.h"

Problem *Problem::factory(std::string solver_type) {
    try {

        if (solver_type == "newmark_general") {
            return new NewmarkGeneral;
        } else if(solver_type == "newmark_testing") {
            return new NewmarkTesting;
        } else {
            throw std::runtime_error("Runtime Error: Problem type " + solver_type + " not supported.");
        }
    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI::COMM_WORLD.Abort(-1);
        return nullptr;
    }
}


