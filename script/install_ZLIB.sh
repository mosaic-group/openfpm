#! /bin/bash

rm zlib-1.3.1.tar.gz
rm -rf zlib-1.3.1
wget http://ppmcore.mpi-cbg.de/upload/zlib-1.3.1.tar.gz -O zlib-1.3.1.tar.gz
if [ $? -ne 0 ]; then
  echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
  exit 1
fi
tar -xf zlib-1.3.1.tar.gz
cd zlib-1.3.1

CC=mpicc CFLAGS=-fPIC  ./configure --prefix=$1/ZLIB
make -j $2
if [ $? -eq 0 ]; then
  make check install
else
  echo -e "\033[91;5;1m ZLIB Installation FAILED \033[0m"
  exit 1
fi
cd ..
