#!/bin/sh

# PETSc
cd
git clone -b maint-3.6 https://bitbucket.org/petsc/petsc petsc-src;
cd petsc-src

# Configure and install.
./configure --download-exodusii --download-netcdf --download-hdf5 --download-chaco \
            --prefix=/home/travis/petsc
make 
make install

# Eigen 3.2
cd
hg clone https://bitbucket.org/eigen/eigen -r 3.2

# Back into salvus directory
cd
cd build/SalvusHub/salvus
