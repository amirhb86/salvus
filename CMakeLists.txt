cmake_minimum_required(VERSION 2.8.7)
project(salvus)

###############################################
# Set all necessary variables here.
# These variables can be set at the command line:
# `cmake ../ -DPETSC_DIR=/opt/petsc`
# or via `ccmake ../`, hitting `c`, editing the
# respective variable, and then `g`
SET(PETSC_DIR /opt/petsc CACHE PATH "PETSC install directory")
SET(EIGEN_INCLUDE /usr/include/eigen3 CACHE PATH "Eigen install directory")
###############################################

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of the builds, options are: Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -std=c++11")
set(CMAKE_CXX_FLAGS_COVERAGE "-fprofile-arcs -ftest-coverage -O0 -g -std=c++11")
MARK_AS_ADVANCED(
        CMAE_CXX_FLAGS_COVERAGE )

find_package(MPI REQUIRED)

link_directories(${PETSC_DIR}/lib)

FILE(GLOB QuadAutoGen src/cxx/Element/HyperCube/Quad/Autogen/*.c)
FILE(GLOB HexAutoGen src/cxx/Element/HyperCube/Hex/Autogen/*.c)
FILE(GLOB TriAutoGen src/cxx/Element/Simplex/Triangle/Autogen/*.c)
FILE(GLOB TetAutoGen src/cxx/Element/Simplex/Tetrahedra/Autogen/*.c)

include_directories(
        ${MPI_INCLUDE_PATH}
        ${PETSC_DIR}/include/
        ${EIGEN_INCLUDE}
        include/)

set(SALVUS_SOURCES
        src/cxx/Mesh/Mesh.cpp
        src/cxx/Problem/Problem.cpp
        src/cxx/Problem/Order2Newmark.cpp
        src/cxx/Element/Simplex/Triangle.cpp
        src/cxx/Element/Simplex/Triangle/TriP1.cpp
        src/cxx/Element/Simplex/Tetrahedra.cpp
        src/cxx/Element/Simplex/Tetrahedra/TetP1.cpp
        src/cxx/Model/ExodusModel.cpp
        src/cxx/Utilities/Utilities.cpp
        src/cxx/Utilities/Options.cpp
        src/cxx/Source/Source.cpp
        src/cxx/Receiver/Receiver.cpp
        src/cxx/Receiver/ReceiverHdf5.cpp
        src/cxx/Element/HyperCube/TensorQuad.cpp
        src/cxx/Element/HyperCube/Quad/QuadP1.cpp
        src/cxx/Element/HyperCube/Hexahedra.cpp
        src/cxx/Element/HyperCube/Hex/HexP1.cpp
        src/cxx/Physics/Element/Scalar.cpp
        src/cxx/Physics/Element/Elastic2D.cpp
        src/cxx/Physics/Element/Elastic3D.cpp
        src/cxx/Physics/Boundary/HomogeneousDirichlet.cpp
        src/cxx/Utilities/kdtree.c
        src/cxx/Utilities/Logging.cpp
        src/cxx/Element/Element.cpp
        src/cxx/Element/HyperCube/Hexahedra.cpp
        src/cxx/Physics/Coupling/AcousticToElastic2D.cpp
        src/cxx/Physics/Coupling/ElasticToAcoustic.cpp
        ${QuadAutoGen}
        ${HexAutoGen}
        ${TriAutoGen}
        ${TetAutoGen}
        )

# build most of our project as a "library", which is shared between
# the main executable and the testing suite.
add_library(salvusCommon ${SALVUS_SOURCES})

add_executable(salvus
        src/cxx/Main.cpp)

add_executable(salvus_test
        src/cxx/Testing/test_main.cpp
        src/cxx/Testing/test_QuadP1.cpp
        src/cxx/Testing/test_TensorQuad.cpp
        src/cxx/Testing/test_source.cpp
        src/cxx/Testing/test_HexP1.cpp
        src/cxx/Testing/test_Hex_Scalar_3D.cpp
        src/cxx/Testing/test_Quad_Scalar_2D.cpp
        src/cxx/Testing/test_TensorHex.cpp
        src/cxx/Testing/test_receiver.cpp
        src/cxx/Testing/test_Element.cpp
        src/cxx/Testing/test_Mesh.cpp
        src/cxx/Testing/test_Model.cpp)

# run `make` (or `make -j5`) to build `salvus` executable
target_link_libraries(salvus salvusCommon ${MPI_LIBRARIES} petsc exodus netcdf hdf5)
# run `make salvus_test` to build testing executable.
target_link_libraries(salvus_test salvusCommon ${MPI_LIBRARIES} petsc exodus netcdf hdf5)
set_target_properties(salvus_test PROPERTIES EXCLUDE_FROM_ALL TRUE) # doesn't get built by `make`
