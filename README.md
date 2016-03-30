# salvus
[![Build Status](https://travis-ci.org/SalvusHub/salvus.svg?branch=master)](https://travis-ci.org/SalvusHub/salvus)

Ask not what your spectral element wave propagator can do for you, but what you can do for your spectral element wave propagator.

## Installation

Download PETSc (at least 3.6.x) from http://www.mcs.anl.gov/petsc/download/, unpack it, and install it with all the required libraries. Adjust the `prefix` to where you want it installed.


```bash
$ ./configure --prefix=/opt/petsc --download-exodusii --download-netcdf --download-hdf5 --download-chaco
```

At the end of each command it tells you to run some other command. Do that
until it is done with everything.

You also need to install version 3.x of the `eigen` library and version 3.x of `cmake`:

```bash
$ brew install eigen cmake          # OSX
```

and on Ubuntu 14.04 (which doesn't include cmake 3.x natively)

``` bash
$ sudo apt-get install libeigen3-dev 
$ sudo add-apt-repository ppa:george-edison55/cmake-3.x
$ sudo apt-get update
$ sudo apt-get install cmake
```

Then copy the `CMakeLists.txt.TEMPLATE` file

```bash
cp CMakeLists.txt.TEMPLATE CMakeLists.txt
```

and edit the following block in `CMakeLists.txt`. `PETSC_DIR` is the `prefix` chosen above.

```cmake
###############################################
# Set all necessary variables here.
# Please only try to edit in this block. If something is
# missing add it.
SET(CMAKE_C_COMPILER /usr/bin/mpicc)
SET(CMAKE_CXX_COMPILER /usr/bin/mpicxx)
SET(PETSC_DIR /opt/petsc)
SET(EIGEN_INCLUDE /usr/include/eigen3)
###############################################
```

It is a good idea to do an "out of source build" to keep the `salvus` directory clean of build artifacts.

``` bash
mkdir build
cd build
cmake ../
```

Finally compile `salvus` with

```bash
$ make -j4
```

### Running it

We aim to collect a number of examples including all required data here: https://github.com/SalvusHub/salvus_data
