#include <Problem/Problem.h>
#include <Utilities/Utilities.h>
#include <Problem/NewmarkGeneral.h>

Problem *Problem::factory(std::string solver_type) {
    try {

        if (solver_type == "newmark_general") {
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


