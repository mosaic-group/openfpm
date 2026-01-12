#! /bin/bash

# check if the directory $1/PETSC exist

CXX=$4
CC=$3
F77=$5
FC=$6
CUDA_SUPPORT=$7

function test_configure_options() {
  cd petsc-3.23.7
  $python_command ./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc  --with-cxx-dialect=C++11 $petsc_openmp --with-mpi-dir=$mpi_dir $configure_options2 --with-debugging=0
  error=$?
  cd ..
}

function haveProg() {
    [ -x "$(command -v $1)" ]
}

python_command=python3

wget http://ppmcore.mpi-cbg.de/upload/petsc-3.23.7.tar.gz -O petsc-3.23.7.tar.gz
tar -xf petsc-3.23.7.tar.gz

touch petsc-3.23.7/config/utils/__init__.py

## If some dependencies has been installed feed them to PETSC

configure_options="--with-64-bit-indices "

if [ "$CUDA_SUPPORT" = true ] ; then
	configure_options="$configure_options --with-cuda --with-cudac=nvcc --with-cuda-arch=$8"
fi

if [ -f "$1/METIS/lib/libmetis.so" ]; then
  configure_options="$configure_options --with-metis-include=$1/METIS/include --with-metis-lib=$1/METIS/lib/libmetis.so "
elif [ -f "$1/METIS/lib/libmetis.dylib" ]; then
  configure_options="$configure_options --with-metis-include=$1/METIS/include --with-metis-lib=$1/METIS/lib/libmetis.dylib "
fi

if [ -f "$1/PARMETIS/lib/libparmetis.so" ]; then
  configure_options="$configure_options --with-parmetis-include=$1/PARMETIS/include --with-parmetis-lib=$1/PARMETIS/lib/libparmetis.so "
elif [ -f "$1/PARMETIS/lib/libparmetis.dylib" ]; then
  configure_options="$configure_options --with-parmetis-include=$1/PARMETIS/include --with-parmetis-lib=$1/PARMETIS/lib/libparmetis.dylib "
fi


if [ -d "$1/BOOST" ]; then

  ### We check incrementaly the options
  configure_options="$configure_options --with-boost=yes --with-boost-dir=$1/BOOST "
fi

if [ -d "$1/MPI" ]; then
  mpi_dir="$1/MPI"
else
  mpi_dir=$(dirname "$(dirname "$(which mpic++)")")
fi

configure_options="$configure_options --with-blas-lib=$1/OPENBLAS/lib/libopenblas.a --with-lapack-lib=$1/OPENBLAS/lib/libopenblas.a --with-suitesparse=yes --with-suitesparse-dir=$1/SUITESPARSE "

configure_options2="$configure_options --download-scalapack --download-mumps --download-superlu_dist --download-hypre "
test_configure_options

if [ $error -eq 0 ]; then
	  echo "SCALAPACK, MUMPS, SUPERLU, HYPRE work with PETSC"
	    configure_options="$configure_options --download-scalapack --download-mumps --download-superlu_dist --download-hypre "
fi

rm -rf petsc-3.23.7
tar -xf petsc-3.23.7.tar.gz
touch petsc-3.23.7/config/utils/__init__.py
cd petsc-3.23.7

if [ x"$CXX" != x"icpc" ]; then

  echo "./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc --with-cxx-dialect=C++11 $petsc_openmp  --with-mpi-dir=$mpi_dir $configure_options  --prefix=$1/PETSC --with-debugging=0"
  if [[ "$OSTYPE" != "cygwin" ]]; then
      $python_command ./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc  --with-cxx-dialect=C++11 $petsc_openmp --with-mpi-dir=$mpi_dir $configure_options --prefix=$1/PETSC --with-debugging=0
  else
    echo "Sorry PETSC installation in not supported on CYGWIN"
  fi

else
  echo "./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc --with-cxx-dialect=C++11 $petsc_openmp  --with-mpi-dir=$mpi_dir $configure_options  --prefix=$1/PETSC --with-debugging=0"
  $python_command ./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc  --with-cxx-dialect=C++11 $petsc_openmp --with-mpi-dir=$mpi_dir $configure_options --prefix=$1/PETSC --with-debugging=0
fi

make all -j $2
make install
