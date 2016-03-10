//
// Created by Michael Afanasiev on 2016-01-29.
//

#ifndef SALVUS_UTILITIES_H
#define SALVUS_UTILITIES_H

#include <iosfwd>
#include <string>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

namespace utilities {

    void print_from_root_mpi(const std::string msg);

    int broadcastInt(int send_buffer);
    std::vector<double> broadcastStdVecFromRoot(std::vector<double> &send_buffer);
    std::string broadcastStringFromRoot(std::string send_str);
    std::vector<std::string> broadcastStringVecFromFroot(std::vector<std::string> &send_buffer);

    class PrintFromRoot {
        
    protected:
        std::ostringstream os;
    private:
        PrintFromRoot(const PrintFromRoot&);
        PrintFromRoot& operator =(const PrintFromRoot&);
        
    public:
        PrintFromRoot();
        virtual ~PrintFromRoot();
        std::ostringstream& Get() { return os; }
    };
    
}

#define PRINT_ROOT() utilities::PrintFromRoot().Get()

enum class AcousticFields {
    displacement = 0, velocity = 1, acceleration = 2, force = 3, acceleration_ = 4, mass_matrix = -1,
    mass_matrix_inverse = -2
};

#endif //SALVUS_UTILITIES_H
