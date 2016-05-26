# salvus
[![Build Status](https://travis-ci.org/SalvusHub/salvus.svg?branch=master)](https://travis-ci.org/SalvusHub/salvus)

Ask not what your spectral element wave propagator can do for you, but
what you can do for your spectral element wave propagator.

## Dependencies

1. You need at least g++ 4.8 (for C++11 features)

    ```bash
    brew install gcc48                   # OSX
    ```
    ``` bash
    sudo apt-get install gcc-4.8 g++-4.8 # Ubuntu 12.04
    ```

2. CMake: <https://cmake.org>

    ```bash
    brew install cmake          # OSX
    ```

    ``` bash
    sudo apt-get install cmake  # Ubuntu 12.04
    ```

3. Eigen: <https://eigen.tuxfamily.org>  
    You will need to install version 3.x of the     `eigen` library (we use 3.2.x for the automated builds and tests)

    ```bash
    brew install eigen                  # OSX
    ```
    ``` bash
    sudo apt-get install libeigen3-dev  # Ubuntu 12.04
    ```

4. MPI (Message Passing Interface):  
    We recommend the MPICH implementation <https://www.mpich.org>

5. PETSc: <https://www.mcs.anl.gov/petsc>  
    You need version 3.6.x. Download it from <http://www.mcs.anl.gov/petsc/download>, unpack it, and install it with all the required libraries. 
    Adjust the `prefix` to where you want it installed.
    Salvus requires the following additional packages be used with PETSc 

      * MPI, ExodusII, HDF5, NetCDF
      * A graph partitioner (chaco, ptscotch, parmetis)

	We recommnend using the chaco partitioner. Note that only ptscotch and parmetis provide support for 64-bit indices.

	PETSc can instal the above packages for you, or use pre-existing installations. See below for additional notes.

    **PETSc installation on a local system without a batch queuing system**  
    If you do not have MPI installed on your machine, PETSc can install it for you
    
    ``` bash
    ./configure --prefix=/opt/petsc \  
      --with-cc=gcc-4.8 --with-cxx=g++-4.8 \  
      --download-mpich=yes --download-exodusii=yes \  
      --download-netcdf=yes --download-hdf5=yes \  
      --download-chaco=yes
    ```

    **PETSc installation on systems with a batch queuing system**  
    You will need to point PETSc to a working MPI implementation provided by target system.
    All other required packages can be installed by PETSc, or you can use local installations 
    provided by your target system if they are available (e.g. HDF5, NetCDF)

    Suppose your target system provides HDF5 and NetCDF (together with MPI), 
    then you would configure PETSc like this

    ``` bash
    ./configure --prefix=/home/software/petsc \  
      --with-batch=no \  
      --with-cc=/path/to/mpicc --with-cxx=/path/to/mpicxx \  
      --with-mpi-dir=/path/to/mpi \  
      --with-netcdf-dir=/path/to/netcdf --with-hdf5-dir=/path/to/h5 \  
      --download-exodusii=yes --download-chaco=yes
    ```

    After a successful configure, follow the instructions issued by PETSc.

    If you have problems configuring PETSc, please refer here <http://www.mcs.anl.gov/petsc/documentation/installation.html>.  
    For serious problems which cannot be resolved, email  
    the entire configure.log and make.log files (as  attachments) to <petsc-maint@mcs.anl.gov>


## Installation

It is a good idea to do an "out of source build" to keep the `salvus`
directory clean of build artifacts.

``` bash
mkdir build
cd build
```

Several variables will need setting in order to successfully compile. From you build directory, you can directly set the `PETSC_DIR`, and `EIGEN_INCLUDE` directories via `ccmake ../`. Hit the `c` key to "configure", make your changes, and hit `g` to generate the Makefiles. Alternatively, you can achieve this via `cmake` on the command line

``` bash
CC=/opt/petsc/bin/mpicc CXX=/opt/petsc/bin/mpicxx cmake ../ -DPETSC_DIR=/opt/petsc -DEIGEN_INCLUDE=/usr/include/eigen3
```

Note the usage of `CC=/opt/petsc/bin/mpicc` `CXX=/opt/petsc/bin/mpicxx`, which is used to change the default compiler used. This only works the **first** time `cmake` is run (it gets cached). 

By default, Salvus is compiled with optimizations (i.e., a release build). To compile for debugging (which adds `-g` and removes `-O3`), add `-DCMAKE_BUILD_TYPE=Debug` (instead of `Release`) to the **first** run of `cmake ../`.

Finally compile `salvus` with

```bash
$ make -j4
```

## Verifying installation

Describe how a test suite can be executed following installation


## Running it

We aim to collect a number of examples including all required data here: <https://github.com/SalvusHub/salvus_data>
