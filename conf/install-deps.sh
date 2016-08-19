#!/bin/sh

# PETSc
cd
# check if petsc has been installed before
if [ -f /home/travis/petsc/include/petsc.h ]; then
    echo "Using cached PETSC install"
else    
    # petsc not in cache, needs installing (or reinstalling)
    echo "Installing petsc to $HOME/petsc"    
    # git clone -b maint https://bitbucket.org/petsc/petsc petsc-src;
    # currently needed for new closure mapping until it is merged into main
    git clone -b maint https://bitbucket.org/petsc/petsc petsc-src
    cd petsc-src

    # Get new branch with fixed spectral closure
    git fetch && git checkout knepley/fix-plex-spectral-order

    # Configure and install.
    ./configure --download-exodusii --download-netcdf --download-hdf5 --download-chaco \
                --prefix=/home/travis/petsc
    make 
    # $HOME/petsc will be cached
    make install
fi
# Eigen 3.2
cd
hg clone https://bitbucket.org/eigen/eigen -r 3.2

# Back into salvus directory
cd
cd build/SalvusHub/salvus
