cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm LANGUAGES C CXX)

if (POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

set(openfpm_VERSION 5.0.0)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake_modules/ /opt/rocm/hip/cmake)

set(BOOST_INCLUDE ${Boost_INCLUDE_DIR} CACHE PATH "Include directory for BOOST")
set(PARMETIS_ROOT CACHE PATH "Parmetis root directory")
set(METIS_ROOT CACHE PATH "Metis root directory")
set(LIBHILBERT_ROOT CACHE PATH "LibHilbert root path")
set(HDF5_ROOT CACHE PATH "HDF5 root path")
set(TINYOBJLOADER_ROOT CACHE PATH "TinyObjLoader library path")
set(SE_CLASS1 CACHE BOOL "Activate compilation with SE_CLASS1")
set(SE_CLASS2 CACHE BOOL "Activate compilation with SE_CLASS2")
set(SE_CLASS3 CACHE BOOL "Activate compilation with SE_CLASS3")
set(ACTION_ON_ERROR CACHE STRING "Action to perform in case of error")
set(PROFILE_WITH_SCOREP CACHE BOOL "Enable profiling with scorep")
set(ALPAKA_ROOT CACHE PATH "Alpaka root path")
set(CPACK_RUN_INSTALL_DEPENDENCIES CACHE BOOL "Set to true if we are creating deb or RPM packages")
set(ENABLE_GARBAGE_INJECTOR CACHE BOOL "Enable the injector of garbage in the memory allocator")
set(ENABLE_VCLUSTER_GARBAGE_INJECTOR CACHE BOOL "Enable the injector of garbage in the vcluster memory buffers")
set(HIP_ENABLE CACHE BOOL "Enable HIP compiler")
set(AMD_ARCH_COMPILE "gfx900" CACHE STRING "AMD gpu architecture used to compile kernels")
set(TEST_PERFORMANCE CACHE BOOL "Enable test performance")
set(CUDA_ON_BACKEND CACHE STRING "Activate several backend for CUDA")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3")

set(ENV{BOOST_ROOT} ${BOOST_ROOT})
set(ENV{HDF5_ROOT} ${HDF5_ROOT})
set(METIS_DIR ${METIS_ROOT})
set(PARMETIS_DIR ${PARMETIS_ROOT})
set(Boost_NO_BOOST_CMAKE OFF)
set(Vc_DIR "${Vc_ROOT}/lib/cmake/Vc/")
message("Searching Vc in ${Vc_DIR}")

find_package(Boost 1.72.0 COMPONENTS unit_test_framework iostreams program_options system filesystem thread OPTIONAL_COMPONENTS fiber context)
find_package(LibHilbert REQUIRED)
find_package(Vc REQUIRED)
find_package(HDF5 REQUIRED)
find_package(TinyObjLoader)
find_package(MPI REQUIRED)
find_package(METIS)
find_package(ParMETIS)

if (CUDA_ON_BACKEND STREQUAL "HIP")
	find_package(HIP)
endif()

if (NOT CUDA_ON_BACKEND STREQUAL "HIP")
	find_package(OpenMP)
endif()

if (CUDA_ON_BACKEND STREQUAL "CUDA")
	find_package(CUDA)
	enable_language(CUDA)
endif()

if (ENABLE_NUMERICS)
	set(ENV{PETSC_DIR} ${PETSC_ROOT})
	set(ENV{EIGEN3_ROOT} ${EIGEN3_ROOT})
	set(SuiteSparse_ROOT ${SUITESPARSE_ROOT})

	set(PETSC_ROOT CACHE PATH "If compiling with linear algebra indicate the PETSC root directory")
	set(EIGEN3_ROOT CACHE PATH "Eigen3 include path")
	set(BLITZ_ROOT CACHE PATH "Blitz root directory")
	set(OPENBLAS_ROOT CACHE PATH "Root path for blas library")
	set(ALGOIM_ROOT CACHE PATH "Algoim root directory")

	find_package(PETSc)
	find_package(Eigen3)
	find_package(BLAS)
	find_package(LAPACK)
	find_package(SuiteSparse OPTIONAL_COMPONENTS UMFPACK)
endif()

if(CUDA_FOUND OR HIP_FOUND)
	set(OPENFPM_INIT_FILE "initialize/initialize_wrapper_cuda.cu")
else()
	set(OPENFPM_INIT_FILE "initialize/initialize_wrapper_cpu.cpp")
endif()

set(CMAKE_SKIP_BUILD_RPATH TRUE)

if(PROFILE_WITH_SCOREP)
	set(CMAKE_CXX_COMPILER_LAUNCHER "scorep")
	set(CMAKE_CC_COMPILER_LAUNCHER "scorep")
	set(CMAKE_CUDA_COMPILER_LAUNCHER "scorep")
endif()

if(ACTION_ON_ERROR)
	set(DEFINE_ACTION_ON_ERROR "#define ${ACTION_ON_ERROR}")
endif()

if(SE_CLASS1)
	set(DEFINE_SE_CLASS1 "#define SE_CLASS1")
endif()

if(SE_CLASS2)
	set(DEFINE_SE_CLASS2 "#define SE_CLASS2")
endif()

if(SE_CLASS3)
	set(DEFINE_SE_CLASS3 "#define SE_CLASS3")
endif()

if(SCAN_COVERTY)
	set(DEFINE_SCAN_COVERTY "#define COVERTY_SCAN")
endif()

if(TEST_PERFORMANCE)
	set(DEFINE_PERFORMANCE_TEST "#define PERFORMANCE_TEST")
endif()

if (TEST_COVERAGE)
	set(DEFINE_TEST_COVERAGE_MODE "#define TEST_COVERAGE_MODE")
endif()

if(APPLE)
	set(DEFINE_HAVE_OSX "#define HAVE_OSX")
endif()

if (Boost_FOUND)
	# ignore BOOST deprecated headers
	add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")

	if ((CUDA_ON_BACKEND STREQUAL "SEQUENTIAL" OR CUDA_ON_BACKEND STREQUAL "OpenMP") AND NOT Boost_CONTEXT_FOUND)
		message( FATAL_ERROR "BOOST is invalid reinstalling" )
	endif()
	set(DEFINE_HAVE_BOOST "#define HAVE_BOOST")
	set(DEFINE_HAVE_BOOST_IOSTREAMS "#define HAVE_BOOST_IOSTREAMS")
	set(DEFINE_HAVE_BOOST_PROGRAM_OPTIONS "#define HAVE_BOOST_PROGRAM_OPTIONS")
	set(DEFINE_HAVE_BOOST_UNIT_TEST_FRAMEWORK "#define HAVE_BOOST_UNIT_TEST_FRAMEWORK")

	if (Boost_CONTEXT_FOUND)
		set(DEFINE_HAVE_BOOST_CONTEXT "#define HAVE_BOOST_CONTEXT")
		set(OPTIONAL_BOOST_LIBS "-lboost_context")
	else()
		#if context is not there CUDA_ON_CPU cannot be activated
		set(CUDA_ON_CPU OFF)
	endif()

	if (Boost_FIBER_FOUND)
		set(DEFINE_HAVE_BOOST_FIBER "#define HAVE_BOOST_FIBER")
		string(CONCAT OPTIONAL_BOOST_LIBS ${OPTIONAL_BOOST_LIBS}  " -lboost_fiber")
	endif()
else()
	message( FATAL_ERROR "BOOST is required in order to install OpenFPM" )
endif()

add_subdirectory(openfpm_devices)
add_subdirectory(openfpm_data)
add_subdirectory(openfpm_vcluster)
add_subdirectory(openfpm_io)
add_subdirectory(openfpm_pdata)

if (ENABLE_NUMERICS)
	add_subdirectory(openfpm_numerics)
	set(DEFINE_ENABLE_NUMERICS "#define ENABLE_NUMERICS")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/config_cmake.h.in ${PROJECT_BINARY_DIR}/config/config.h)

# WRITE BUILD OPTIONS TO CONFIG FILES. LATER USED TO CREATE example.mk
if (OPENMP_FOUND)
	if (NOT OpenMP_CXX_INCLUDE_DIRS STREQUAL "")
		file(WRITE openmp_flags "${OpenMP_CXX_FLAGS} -I${OpenMP_CXX_INCLUDE_DIRS}")
	else()
		file(WRITE openmp_flags "${OpenMP_CXX_FLAGS}")
	endif()
	string(REPLACE ";" " " OpenMP_CXX_LIBRARIES_LIB "${OpenMP_CXX_LIBRARIES}")
	file(WRITE openmp_libs ${OpenMP_CXX_LIBRARIES_LIB})
else()
	file(WRITE openmp_flags " ")
	file(WRITE openmp_libs " ")
endif()

file(WRITE cuda_lib "${CUDA_cudart_static_LIBRARY} ${CUDA_cudadevrt_LIBRARY}") 
file(WRITE cuda_options "${WARNING_SUPPRESSION_AND_OPTION_NVCC_TEXT}")

file(WRITE mpi_include "-I${MPI_C_INCLUDE_DIRS}")
file(WRITE mpi_libs "${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES}")

if(CUDA_ON_BACKEND STREQUAL "SEQUENTIAL" OR CUDA_ON_BACKEND STREQUAL "OpenMP")
	file(WRITE cuda_include "-D__NVCC__ -DCUDART_VERSION=11000")
	file(WRITE cuda_on_cpu "YES")
else()
	if (NOT CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES STREQUAL "")
		file(WRITE cuda_include "-I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
	else()
		file(WRITE cuda_include " ")
	endif()
	file(WRITE cuda_on_cpu "NO")
endif()

if (Boost_FOUND)
	if (Boost_CONTEXT_FOUND)
	       file(WRITE optional_boost_libs "-L${BOOST_ROOT}/lib -lboost_context")
	else()
	       file(WRITE optional_boost_libs " ")
	endif()
endif()

if(HIP_FOUND)
	file(WRITE hip_enabled 1)
else()
	file(WRITE hip_enabled 0)
endif()

#################### CPack to create auto installing packages

include(InstallRequiredSystemLibraries)

string(REPLACE "." ";" VERSION_LIST ${openfpm_VERSION})
list(GET VERSION_LIST 0 OPENFPM_VERSION_MAJOR)
list(GET VERSION_LIST 1 OPENFPM_VERSION_MINOR)
list(GET VERSION_LIST 2 OPENFPM_VERSION_PATCH)

if (CPACK_RUN_INSTALL_DEPENDENCIES)

	###### Fix post inst script ######

	

	######                      ######

	set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenFPM distributed data-structures")
	set(CPACK_PACKAGE_VENDOR "Vector Cluster")
	set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
	set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
	set(CPACK_PACKAGE_VERSION_MAJOR ${OPENFPM_VERSION_MAJOR})
	set(CPACK_PACKAGE_VERSION_MINOR ${OPENFPM_VERSION_MINOR})
	set(CPACK_PACKAGE_VERSION_PATCH ${OPENFPM_VERSION_PATCH})
	set(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
	set(CPACK_PACKAGE_INSTALL_DIRECTORY /usr/local/openfpm)
	set(CPACK_PACKAGING_INSTALL_PREFIX /usr/local/openfpm)
	set(CPACK_PACKAGE_HOMEPAGE_URL http://openfpm.mpi-cbg.de)
	set(CPACK_RPM_PACKAGE_AUTOREQPROV NO)
	set(CPACK_PACKAGE_CONTACT abhinavsns7@gmail.com)
	set(CPACK_DEBIAN_PACKAGE_MAINTAINER Pietro Incardona)
	set(CPACK_PKG_PACKAGE_MAINTAINER Abhinav Singh)
	set(CPACK_RPM_POST_INSTALL_SCRIPT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/src/scripts/postinst)
	set(CPACK_POSTFLIGHT_OPENFPM_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/src/scripts/postflight)
	set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA ${CMAKE_CURRENT_SOURCE_DIR}/src/scripts/postinst)
	set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")

	install(FILES $ENV{DEP_PACKING}/openfpm_vars
			      $ENV{DEP_PACKING}/.openfpm.mk
		DESTINATION source
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/BOOST
		DESTINATION dependencies/
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/EIGEN
		DESTINATION dependencies/
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/HDF5
		DESTINATION dependencies/
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/LIBHILBERT
		DESTINATION dependencies/
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/METIS
		DESTINATION dependencies/
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/MPI
		DESTINATION dependencies/
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/OPENBLAS
		DESTINATION dependencies/
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/PARMETIS
		DESTINATION dependencies/
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/PETSC
		DESTINATION dependencies/
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/SUITESPARSE
		DESTINATION dependencies/
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/VCDEVEL
		DESTINATION dependencies/
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/BLITZ
		DESTINATION dependencies/
		COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/ALGOIM
		DESTINATION dependencies/
		COMPONENT OpenFPM)

endif()

include(CPack)

cpack_add_component(OpenFPM
		    DISPLAY_NAME OpenFPM
		    DESCRITION OpenFPM binary files)
