#! /bin/bash

# check if the directory $1/PETSC exist

CXX=$4
CC=$3
F77=$5
FC=$6

function test_configure_options() {
  cd petsc-3.20.5
  $python_command ./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc  --with-cxx-dialect=C++11 $petsc_openmp --with-mpi-dir=$mpi_dir $configure_options2 --with-debugging=0
  error=$?
  cd ..
}

function haveProg() {
    [ -x "$(command -v $1)" ]
}

python_command=python3

wget http://ppmcore.mpi-cbg.de/upload/petsc-3.20.5.tar.gz -O petsc-3.20.5.tar.gz
tar -xf petsc-3.20.5.tar.gz

## If some dependencies has been installed feed them to PETSC

configure_options="--with-64-bit-indices --with-parmetis-include=$1/PARMETIS/include --with-parmetis-lib=$1/PARMETIS/lib/libparmetis.a"

if [ -f "$1/METIS/lib/libmetis.so" ]; then
  configure_options="$configure_options --with-metis-include=$1/METIS/include --with-metis-lib=$1/METIS/lib/libmetis.so "
elif [ -f "$1/METIS/lib/libmetis.dylib" ]; then
  configure_options="$configure_options --with-metis-include=$1/METIS/include --with-metis-lib=$1/METIS/lib/libmetis.dylib "
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


configure_options="$configure_options --with-blas-lib=$1/OPENBLAS/lib/libopenblas.a --with-lapack-lib=$1/OPENBLAS/lib/libopenblas.a"

echo "Testing if PETSC work with SUITESPARSE"
configure_options2="$configure_options --with-suitesparse=yes --with-suitesparse-dir=$1/SUITESPARSE "
test_configure_options

if [ $error -eq 0 ]; then
  echo "SUITESPARSE work with PETSC"
  configure_options="$configure_options --with-suitesparse=yes --with-suitesparse-dir=$1/SUITESPARSE "
fi

configure_options2="$configure_options --download-scalapack"
test_configure_options

if [ $error -eq 0 ]; then
  echo "SCALAPACK work with PETSC"
  configure_options="$configure_options --download-scalapack "
fi

configure_options2="$configure_options --download-mumps"
test_configure_options

if [ $error -eq 0 ]; then
  echo "MUMPS work with PETSC"
  configure_options="$configure_options --download-mumps"
fi

echo "Testing if PETSC work with SUPERLU"
configure_options2="$configure_options --download-superlu_dist "
test_configure_options

if [ $error -eq 0 ]; then
  echo "SUPERLU work with PETSC"
  configure_options="$configure_options --download-superlu_dist "
fi

configure_options2="$configure_options --download-hypre"
test_configure_options

if [ $error -eq 0 ]; then
  echo "HYPRE work with PETSC"
  configure_options="$configure_options --download-hypre"
fi

configure_options="$configure_options --download-scalapack "

rm -rf petsc-3.20.5
tar -xf petsc-3.20.5.tar.gz
cd petsc-3.20.5

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
