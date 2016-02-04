cmake_minimum_required(VERSION 3.2.2)
project(salvus)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
SET(CMAKE_C_COMPILER /usr/local/bin/mpicc)
SET(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules")

SET(PETSC_DIR PATH)
find_package(PETSc REQUIRED CXX)
include_directories("${PETSC_DIR}/include")
link_directories(${PETSC_DIR}/lib)

include_directories(
        /usr/local/include/eigen3
        cxx
        cxx/Element
        cxx/Element/HyperCube
        cxx/Element/HyperCube/Square
        cxx/Element/HyperCube/Square/Autogen
        cxx/Element/HyperCube/Cube
        cxx/Element/Simplex
        cxx/Element/Simplex/Triangle
        cxx/Element/Simplex/Tetrahedron
        cxx/Model)

add_executable(salvus
        src/cxx/Main.cpp
        src/cxx/Solver.cpp
        src/cxx/Mesh.cpp
        src/cxx/Utilities.cpp
        src/cxx/Options.cpp
        src/cxx/Element/Element.cpp
        src/cxx/Element/HyperCube/Square.cpp
        src/cxx/Element/HyperCube/Square/SquareAcousticOrderFour.cpp
        src/cxx/Element/HyperCube/Square/Autogen/order4_square.c
        src/cxx/Model/ExodusModel.cpp
        src/cxx/kdtree.c src/cxx/Source.cpp src/cxx/Source.h)

target_link_libraries(salvus petsc exodus netcdf)
