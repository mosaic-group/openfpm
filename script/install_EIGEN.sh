#! /bin/bash

wget http://ppmcore.mpi-cbg.de/upload/eigen-3.4.0.tar.bz2 -O eigen-3.4.0.tar.bz2
tar -xf eigen-3.4.0.tar.bz2

cd eigen-3.4.0/
rm -rf $1/EIGEN/
mkdir $1/EIGEN/
mv Eigen $1/EIGEN/Eigen
touch $1/EIGEN/signature_of_eigen3_matrix_library
