#! /bin/bash

rm -rf blitz-1.0.2
wget http://ppmcore.mpi-cbg.de/upload/blitz-1.0.2.tar.gz -O blitz-1.0.2.tar.gz
tar -xf blitz-1.0.2.tar.gz
cd blitz-1.0.2

BUILDDIR=build
mkdir -p $BUILDDIR
cd $BUILDDIR
echo "cmake ../. -DCMAKE_INSTALL_PREFIX=$1/BLITZ"
cmake ../. -DCMAKE_INSTALL_PREFIX=$1/BLITZ -DBZ_HAVE_LIBPAPI=OFF
make -j $2
make install
