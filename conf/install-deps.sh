#!/bin/sh

# PETSc
cd
git clone -b maint https://bitbucket.org/petsc/petsc petsc-src;
cd petsc-src

# Configure and install.
./configure --download-exodusii --download-netcdf --download-hdf5 --download-chaco \
            --prefix=/home/travis/petsc
make 
make install

# Eigen
cd
hg clone https://bitbucket.org/eigen/eigen

# Back into salvus directory
cd
cd build/SalvusHub/salvus
