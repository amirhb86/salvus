#!/bin/sh

# PETSc
cd
git clone -b maint https://bitbucket.org/petsc/petsc petsc;
cd petsc

# Configure and install.
./configure --download-exodusii --download-netcdf --download-hdf5 --download-chaco \
--with-cc=/usr/bin/gcc-4.8 --with-fc=/usr/bin/gfortran-4.8 --with-cxx=/usr/bin/g++-4.8 \
--download-openmpi=yes
make PETSC_DIR=/home/travis/petsc PETSC_ARCH=arch-linux2-c-debug all
make PETSC_DIR=/home/travis/petsc PETSC_ARCH=arch-linux2-c-debug test
make PETSC_DIR=/home/travis/petsc PETSC_ARCH=arch-linux2-c-debug streams NPMAX=2

# Eigen
cd
hg clone https://bitbucket.org/eigen/eigen

# Back into salvus directory
cd
cd build/SalvusHub/salvus
