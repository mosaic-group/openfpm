#! /bin/bash

rm -rf SuiteSparse
wget http://ppmcore.mpi-cbg.de/upload/SuiteSparse-7.11.0.tar.gz -O SuiteSparse-7.11.0.tar.gz
tar -xf SuiteSparse-7.11.0.tar.gz
cd SuiteSparse-7.11.0
make clean

if [ x"$CXX" == x"icpc" ]; then
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/$1/OPENBLAS/lib"
    STS_LIB="-shared-intel -lrt -lifcore"
fi

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$1/OPENBLAS/lib"
if [[ "$OSTYPE" == "cygwin" ]]; then
    export PATH="$PATH:$(pwd)/lib"
    echo "$PATH"
fi

export CMAKE_INSTALL_PREFIX="$1/SUITESPARSE"

CMAKE_OPTIONS="-DCMAKE_C_COMPILER=$3 -DCMAKE_CXX_COMPILER=$4 -DBLA_VENDOR=OpenBLAS -DCMAKE_INSTALL_PREFIX=$1/SUITESPARSE -DSUITESPARSE_USE_CUDA=OFF -DCHOLMOD_USE_CUDA=OFF -DSPQR_USE_CUDA=OFF -DBLAS_ROOT=$1/OPENBLAS/ -DLAPACK_ROOT=$1/OPENBLAS/" make VERBOSE=1 library
make install