#! /bin/bash

### 1.8.19 does not compile on CYGWIN
rm -rf hdf5-1.14.5
wget http://ppmcore.mpi-cbg.de/upload/hdf5-1.14.5.tar.gz -O hdf5-1.14.5.tar.gz
tar -xf hdf5-1.14.3.tar.gz
cd hdf5-1.14.3

# Disable zlib is completly unstable
if [[ "$OSTYPE" != "cygwin" ]]; then
        CC=$3 ./configure --with-zlib=$1/ZLIB --enable-parallel --prefix=$1/HDF5
	make -j $2
else
        CC=$3 ./configure --enable-parallel --prefix=$1/HDF5
	make CFLAGS=-D_POSIX_C_SOURCE -j $2
fi

mkdir $1/HDF5
make install
