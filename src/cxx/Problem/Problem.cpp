//
// Created by Michael Afanasiev on 2016-01-27.
//

#include "Problem.h"
#include "TimeDomainElastic2d.h"
#include "TimeDomainScalar2d.h"
#include "NewmarkGeneral.h"

Problem *Problem::factory(std::string solver_type) {
    try {

        if (solver_type == "solver_newmark_dimension_2d_physics_acoustic") {
            return new TimeDomainScalar2d;
        } else if (solver_type == "solver_newmark_dimension_2d_physics_elastic") {
            return new TimeDomainElastic2d;
        } else if (solver_type == "newmark_general") {
            return new NewmarkGeneral;
        } else {
            throw std::runtime_error("Runtime Error: Problem type " + solver_type + " not supported.");
        }
    } catch (std::exception &e) {
        utilities::print_from_root_mpi(e.what());
        MPI::COMM_WORLD.Abort(-1);
        return nullptr;
    }
}


