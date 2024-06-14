#! /bin/bash

# None Debug Release RelWithDebInfo
if [[ -v BUILD_TYPE ]]; then
    build_type="$BUILD_TYPE"
else
    build_type=Release
fi

if [[ -v SE_CLASS1 ]]; then
    se_class1="$SE_CLASS1"
else
    se_class1=OFF
fi

if [[ -v SE_CLASS2 ]]; then
    se_class2="$SE_CLASS2"
else
    se_class2=OFF
fi

if [[ -v SE_CLASS3 ]]; then
    se_class3="$SE_CLASS3"
else
    se_class3=OFF
fi

if [[ -v TEST_COVERAGE ]]; then
    test_coverage="$TEST_COVERAGE"
else
    test_coverage=OFF
fi

if [[ -v SCAN_COVERTY ]]; then
    scan_coverty="$SCAN_COVERTY"
else
    scan_coverty=OFF
fi

# NONE, CUDA, SEQUENTIAL, OpenMP, HIP
if [[ -v CUDA_ON_BACKEND ]]; then
    cuda_on_backend="$CUDA_ON_BACKEND"
else
    cuda_on_backend=NONE
fi

if [[ -v ENABLE_NUMERICS ]]; then
    enable_numerics="$ENABLE_NUMERICS"
else
    enable_numerics=ON
fi

if [[ -v TEST_PERFORMANCE ]]; then
    test_performance="$TEST_PERFORMANCE"
else
    test_performance=OFF
fi

if [[ -v ENABLE_ASAN ]]; then
    enable_asan="$ENABLE_ASAN"
else
    enable_asan=OFF
fi

if [[ -v ENABLE_GARBAGE_INJECTOR ]]; then
    enable_garbage_injector="$ENABLE_GARBAGE_INJECTOR"
else
    enable_garbage_injector=OFF
fi

if [[ -v ENABLE_VCLUSTER_GARBAGE_INJECTOR ]]; then
    enable_vcluster_garbage_injector="$ENABLE_VCLUSTER_GARBAGE_INJECTOR"
else
    enable_vcluster_garbage_injector=OFF
fi

# -----------

prefix_dependencies="$1"
configure_options="-DCMAKE_INSTALL_PREFIX=$2 "
ld_lib_pathopt=

configure_options="$configure_options -DCMAKE_BUILD_TYPE=$cmake_build_type "
configure_options="$configure_options -DSE_CLASS1=$se_class1 "
configure_options="$configure_options -DSE_CLASS2=$se_class2 "
configure_options="$configure_options -DSE_CLASS3=$se_class3 "
configure_options="$configure_options -DTEST_COVERAGE=$test_coverage "
configure_options="$configure_options -DSCAN_COVERTY=$scan_coverty "
configure_options="$configure_options -DTEST_PERFORMANCE=$test_performance "
configure_options="$configure_options -DENABLE_ASAN=$enable_asan "
configure_options="$configure_options -DENABLE_NUMERICS=$enable_numerics "
configure_options="$configure_options -DENABLE_GARBAGE_INJECTOR=$enable_garbage_injector "
configure_options="$configure_options -DENABLE_VCLUSTER_GARBAGE_INJECTOR=$enable_vcluster_garbage_injector "


configure_options="$configure_options -DCUDA_ON_BACKEND=$cuda_on_backend "

if [ x"$cuda_on_backend" != x"CUDA" ]; then
    configure_options="$configure_options"
else
    configure_options="$configure_options -DCMAKE_CUDA_HOST_COMPILER=$(which $CXX) "
fi
if [ x"$CXXCUDA" == x"" ]; then
    configure_options="$configure_options"
else
    configure_options="$configure_options -DCMAKE_CUDA_COMPILER=$(which $CXXCUDA) "
fi

if [ -d "$prefix_dependencies/MPI" ]; then
    configure_options="$configure_options -DMPI_VENDOR=openmpi "
    configure_options="$configure_options -DMPI_ROOT=$prefix_dependencies/MPI "
fi

if [ -d "$prefix_dependencies/PETSC" ]; then
    configure_options="$configure_options -DPETSC_ROOT=$prefix_dependencies/PETSC "
fi

if [ -d "$prefix_dependencies/BOOST" ]; then
    configure_options=" $configure_options -DBOOST_ROOT=$prefix_dependencies/BOOST -DBoost_NO_BOOST_CMAKE=ON "
fi

if [ -d "$prefix_dependencies/HDF5" ]; then
    configure_options=" $configure_options -DHDF5_ROOT=$prefix_dependencies/HDF5/ "
fi

if [ -d "$prefix_dependencies/LIBHILBERT" ]; then
    configure_options="$configure_options -DLIBHILBERT_ROOT=$prefix_dependencies/LIBHILBERT "
fi

if [ -d "$prefix_dependencies/BLITZ" ]; then
    configure_options=" $configure_options -DBLITZ_ROOT=$prefix_dependencies/BLITZ "
fi

if [ -d "$prefix_dependencies/ALGOIM" ]; then
    configure_options=" $configure_options -DALGOIM_ROOT=$prefix_dependencies/ALGOIM "
fi

if [ -d "$prefix_dependencies/PARMETIS" ]; then
    configure_options=" $configure_options -DPARMETIS_ROOT=$prefix_dependencies/PARMETIS "
fi

if [ -d "$prefix_dependencies/METIS" ]; then
	configure_options=" $configure_options -DMETIS_ROOT=$prefix_dependencies/METIS "
fi

if [ -d "$prefix_dependencies/VCDEVEL" ]; then
    configure_options=" $configure_options -DVc_ROOT=$prefix_dependencies/VCDEVEL "
fi

if [ -d "$prefix_dependencies/OPENBLAS" ]; then
    configure_options=" $configure_options -DOPENBLAS_ROOT=$prefix_dependencies/OPENBLAS/"
fi

if [ -d "$prefix_dependencies/SUITESPARSE"  -a -f "$prefix_dependencies/SUITESPARSE/include/umfpack.h" ]; then
    configure_options="$configure_options -DSUITESPARSE_ROOT=$prefix_dependencies/SUITESPARSE "
    ld_lib_pathopt=$prefix_dependencies/SUITESPARSE/lib
fi

if [ -d "$prefix_dependencies/EIGEN" ]; then
    configure_options=" $configure_options -DEIGEN3_ROOT=$prefix_dependencies/EIGEN "
fi

if [ -d "$i_dir/CATALYST/lib" ]; then
  configure_options=" $configure_options -Dcatalyst_DIR=$prefix_dependencies/CATALYST/lib/cmake/catalyst-2.0 "
elif [ -d "$i_dir/CATALYST/lib64" ]; then
  configure_options=" $configure_options -Dcatalyst_DIR=$prefix_dependencies/CATALYST/lib64/cmake/catalyst-2.0 "
fi

echo "CXX=mpic++ DYLD_LIBRARY_PATH=$ld_lib_pathopt cmake ../. $configure_options"
printf "CXX=mpic++ DYLD_LIBRARY_PATH=$ld_lib_pathopt cmake ../. $configure_options" > cmake_build_options