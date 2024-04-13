#! /bin/bash

prefix_dependencies="$1"
prefix_openfpm="$2"

openmp_libs="$(cat openmp_libs)"
openmp_flags="$(cat openmp_flags)"
cuda_include_dirs="$(cat cuda_include)"
cuda_on_cpu="$(cat cuda_on_cpu)"
cuda_lib="$(cat cuda_lib)"
cuda_options="$(cat cuda_options)"

hip_enabled=$(cat hip_enabled)
if [ x"$hip_enabled" == x"1" ]; then
	mpi_include_dirs=$(cat mpi_include)
	mpi_libs=$(cat mpi_libs)
fi

optional_boost=$(cat optional_boost_libs)

if [ -d "$prefix_dependencies/HDF5/lib" ]; then
  hdf5_lib_dir=-L$prefix_dependencies/HDF5/lib
elif [ -d "$prefix_dependencies/HDF5/lib64" ]; then
  hdf5_lib_dir=-L$prefix_dependencies/HDF5/lib64
fi

lin_alg_dir=""
lin_alg_lib=""
lin_alg_inc=""

if [ -d "$prefix_dependencies/PETSC" ]; then
    lin_alg_dir="$lin_alg_dir -L$prefix_dependencies/PETSC/lib"
    lin_alg_lib="$lin_alg_lib -lpetsc"
    lin_alg_inc="$lin_alg_inc -I$prefix_dependencies/PETSC/include"
fi

if [ -d "$prefix_dependencies/OPENBLAS" ]; then
    lin_alg_dir="$lin_alg_dir -L$prefix_dependencies/OPENBLAS/lib"
    lin_alg_lib="$lin_alg_lib -lopenblas"
    lin_alg_inc="$lin_alg_inc -I$prefix_dependencies/OPENBLAS/include"
fi

if [ -d "$prefix_dependencies/SUITESPARSE"  -a -f "$prefix_dependencies/SUITESPARSE/include/umfpack.h" ]; then
    lin_alg_dir="$lin_alg_dir -L$prefix_dependencies/SUITESPARSE/lib"
    lin_alg_inc="$lin_alg_inc -I$prefix_dependencies/SUITESPARSE/include"
    lin_alg_lib="$lin_alg_lib -lumfpack -lamd -lbtf -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lklu -ldl -lrbio -lspqr -lsuitesparseconfig"
fi

if [ -d "$prefix_dependencies/EIGEN" ]; then
    lin_alg_inc="$lin_alg_inc -I$prefix_dependencies/EIGEN"
fi

if [[ "$OSTYPE" == "linux-gnu" || "$OSTYPE" == "linux" ]]; then
    lin_alg_lib="$lin_alg_lib -lrt"
fi

if [ -d "$prefix_dependencies/CATALYST/lib" ]; then
  catalyst_lib=$prefix_dependencies/CATALYST/lib
  catalyst_lib_dir=-L$prefix_dependencies/CATALYST/lib
elif [ -d "$prefix_dependencies/CATALYST/lib64" ]; then
  catalyst_lib=$prefix_dependencies/CATALYST/lib64
  catalyst_lib_dir=-L$prefix_dependencies/CATALYST/lib64
fi

echo "INCLUDE_PATH=$mpi_include_dirs $cuda_include_dirs $openmp_flags  -I.  -I$prefix_openfpm/openfpm_numerics/include -I$prefix_openfpm/openfpm_pdata/include/config -I$prefix_openfpm/openfpm_pdata/include -I$prefix_openfpm/openfpm_data/include -I$prefix_openfpm/openfpm_vcluster/include -I$prefix_openfpm/openfpm_io/include -I$prefix_openfpm/openfpm_devices/include -I$prefix_dependencies/VCDEVEL/include  -I$prefix_dependencies/METIS/include -I$prefix_dependencies/PARMETIS/include -I$prefix_dependencies/BOOST/include -I$prefix_dependencies/HDF5/include -I$prefix_dependencies/LIBHILBERT/include  $lin_alg_inc -I$prefix_dependencies/BLITZ/include -I$prefix_dependencies/ALGOIM/include  -I$prefix_dependencies/SUITESPARSE/include -I$prefix_dependencies/CATALYST/include/catalyst-2.0 " > example.mk
echo "LIBS_PATH=$mpi_libs -L$prefix_openfpm/openfpm_devices/lib -L$prefix_openfpm/openfpm_pdata/lib  -L$prefix_openfpm/openfpm_vcluster/lib -L$prefix_dependencies/VCDEVEL/lib  -L$prefix_dependencies/METIS/lib -L$prefix_dependencies/PARMETIS/lib  -L$prefix_dependencies/BOOST/lib $hdf5_lib_dir -L$prefix_dependencies/LIBHILBERT/lib  $lin_alg_dir $catalyst_lib_dir " >> example.mk
if [ x"$cuda_on_cpu" == x"YES" ]; then
   echo "CUDA_ON_CPU=YES" >> example.mk
fi
echo "LIBS=$openmp_flags $mpi_libs $openmp_libs -lvcluster -lofpm_pdata -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lboost_program_options -lhdf5 -llibhilbert -lVc  $cuda_lib $lin_alg_lib -ldl -lboost_filesystem -lboost_system $optional_boost -lcatalyst" >> example.mk
echo "LIBS_NVCC=-Xcompiler=$openmp_flags $mpi_libs -lvcluster -lofpm_pdata -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lboost_program_options -lhdf5 -llibhilbert -lVc  $cuda_lib $lin_alg_lib -ldl -lboost_filesystem -lboost_system $optional_boost -lcatalyst" >> example.mk
echo "INCLUDE_PATH_NVCC=-Xcompiler=$openmp_flags "$cuda_options" $mpi_include_dirs -I. -I$prefix_openfpm/openfpm_numerics/include -I$prefix_openfpm/openfpm_pdata/include/config -I$prefix_openfpm/openfpm_pdata/include -I$prefix_openfpm/openfpm_data/include -I$prefix_openfpm/openfpm_vcluster/include -I$prefix_openfpm/openfpm_io/include -I$prefix_openfpm/openfpm_devices/include -I$prefix_dependencies/METIS/include -I$prefix_dependencies/PARMETIS/include -I$prefix_dependencies/BOOST/include -I$prefix_dependencies/HDF5/include -I$prefix_dependencies/LIBHILBERT/include  $lin_alg_inc -I$prefix_dependencies/BLITZ/include -I$prefix_dependencies/ALGOIM/include  -I$prefix_dependencies/SUITESPARSE/include -I$prefix_dependencies/CATALYST/include/catalyst-2.0 " >> example.mk