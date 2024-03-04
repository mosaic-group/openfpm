#! /bin/bash

wget http://ppmcore.mpi-cbg.de/upload/OpenBLAS-0.3.26.tar.gz -O OpenBLAS-0.3.26.tar.gz
tar -xf OpenBLAS-0.3.26.tar.gz
cd OpenBLAS-0.3.26

make FC=$FC CC=$CC -j $2
mkdir $1/OPENBLAS
make install PREFIX=$1/OPENBLAS
