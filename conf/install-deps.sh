#!/bin/sh

# PETSc
cd
# check if petsc has been installed before
if [ ! -f /home/travis/petsc/include/petsc.h ]; then
    echo "Installing petsc to $HOME/petsc"    
    git clone -b maint https://bitbucket.org/petsc/petsc petsc-src;
    cd petsc-src

    # Get Matt's new branch.
    git fetch && git checkout knepley/feature-plex-3d-sem-order

    # Configure and install.
    ./configure --download-exodusii --download-netcdf --download-hdf5 --download-chaco \
                --prefix=/home/travis/petsc
    make 
    # $HOME/petsc will be cached
    make install
else
    echo "Using cached PETSC install"
fi
# Eigen 3.2
cd
hg clone https://bitbucket.org/eigen/eigen -r 3.2

# Back into salvus directory
cd
cd build/SalvusHub/salvus
