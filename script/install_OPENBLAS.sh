#! /bin/bash

wget http://ppmcore.mpi-cbg.de/upload/OpenBLAS-develop.tar.gz -O OpenBLAS-develop.tar.gz
tar -xf OpenBLAS-develop.tar.gz
cd OpenBLAS

make clean
mkdir $1/OPENBLAS
make FC=$6 CC=$3 BUILD_WITHOUT_LAPACK=1 USE_OPENMP=0 NOFORTRAN=1 -j $2
make install PREFIX=$1/OPENBLAS
