#!/bin/sh

# Setup directory tree.
mkdir -p petsc_install
mkdir -p petsc_install/lib
mkdir -p petsc_install/include
mkdir -p petsc_install/bin

# # Linux
# sudo apt-get update -q
#
# # system libs.
# sudo apt-get install build-essential gfortran perl g++ gcc cmake m4 git
# sudo apt-get install -y liblapack-dev liblapack-doc liblapack3gf
# sudo apt-get install -u libblas-dev libblas-doc libblas3gf
# sudo apt-get install -y gfortran openmpi-bin openmpi-common libopenmpi-dev

# PETSc
git clone -b maint https://bitbucket.org/petsc/petsc petsc; cd petsc
./configure --download-exodusii --download-netcdf --download-hdf5 --download-chaco
make PETSC_DIR=/home/travis/petsc PETSC_ARCH=arch-linux2-c-debug all
make PETSC_DIR=/home/travis/petsc PETSC_ARCH=arch-linux2-c-debug test
make PETSC_DIR=/home/travis/petsc PETSC_ARCH=arch-linux2-c-debug streams NPMAX=2
