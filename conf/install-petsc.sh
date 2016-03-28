#!/bin/sh

# PETSc
cd
git clone -b maint https://bitbucket.org/petsc/petsc petsc;
cd petsc

./configure --download-exodusii --download-netcdf --download-hdf5 --download-chaco
make PETSC_DIR=/home/travis/petsc PETSC_ARCH=arch-linux2-c-debug all
make PETSC_DIR=/home/travis/petsc PETSC_ARCH=arch-linux2-c-debug test
make PETSC_DIR=/home/travis/petsc PETSC_ARCH=arch-linux2-c-debug streams NPMAX=2
