cmake_minimum_required(VERSION 3.0.0)
project(salvus)

###############################################
# Set all necessary variables here.
# Please only try to edit in this block. If something is
# missing add it.
SET(CMAKE_C_COMPILER /usr/bin/mpicc)
SET(CMAKE_CXX_COMPILER /usr/bin/mpicxx)
SET(PETSC_DIR /opt/petsc)
SET(EIGEN_INCLUDE /usr/include/eigen3)
###############################################

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

link_directories(${PETSC_DIR}/lib)

FILE(GLOB QuadAutoGen src/cxx/Element/HyperCube/Square/Autogen/*.c)

include_directories(
        ${PETSC_DIR}/include/
        ${EIGEN_INCLUDE}
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
        src/cxx/Element/HyperCube/Square.cpp
        src/cxx/Element/HyperCube/Square/Acoustic.cpp
        src/cxx/Model/ExodusModel.cpp
        src/cxx/kdtree.c
        src/cxx/Source.cpp
        ${QuadAutoGen}
        )

target_link_libraries(salvus petsc exodus netcdf)
