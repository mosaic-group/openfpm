#!/bin/bash

rm -rf openmpi-5.0.8
wget http://ppmcore.mpi-cbg.de/upload/openmpi-5.0.8.tar.gz -O openmpi-5.0.8.tar.gz
tar -xvf openmpi-5.0.8.tar.gz
cd openmpi-5.0.8

if [ x"$3" == x"1" ]; then
   echo "Installing MPI with GPU support"

   # Detect where is nvcc
   cuda_location=$(dirname $(dirname $(which nvcc)) )

   ./configure --with-hwloc=internal --with-libevent=internal $mpi_options --with-cuda=$cuda_location --prefix=$1/MPI --enable-mpi-fortran=yes CC=$4 CXX=$5 F77=$6 FC=$7 $8
else
   echo "Installing MPI without GPU support"
   ./configure --with-hwloc=internal --with-libevent=internal $mpi_options --prefix=$1/MPI --enable-mpi-fortran=yes CC=$4 CXX=$5 F77=$6 FC=$7 $8
fi
make -j $2
make install

