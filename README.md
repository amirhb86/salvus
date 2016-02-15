# salvus
Ask not what your spectral element wave propagator can do for you, but what you can do for your spectral element wave propagator.

## Installation

Download PETSc from http://www.mcs.anl.gov/petsc/download/, unpack it, and install it with all the required libraries. Adjust the `prefix` to where you want it installed.


```bash
$ ./configure --prefix=/opt/petsc --download-exodusii --download-netcdf --download-hdf5 --download-chaco
```

At the end of each command it tells you to run some other command. Do that
until it is done with everything.

You also need to install version 3 of the `eigen` library and `cmake`:

```bash
$ brew install eigen cmake          # OSX
$ sudo apt-get install libeigen3-dev cmake  # Ubuntu
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

Finally compile `salvus` with

```bash
$ cmake .
$ make -j4
```

### Running it

We aim to collect a number of examples including all required data here: https://github.com/SalvusHub/salvus_data