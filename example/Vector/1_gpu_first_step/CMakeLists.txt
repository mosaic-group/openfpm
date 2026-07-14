cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

if (NOT CUDA_ON_BACKEND STREQUAL "METAL")
	message(FATAL_ERROR "The in-tree CMake target is enabled by the Metal backend")
endif()

execute_process(COMMAND "${MPI_CXX_COMPILER}" --showme:incdirs
	OUTPUT_VARIABLE mpi_wrapper_include_text
	OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
separate_arguments(mpi_wrapper_include_dirs NATIVE_COMMAND
	"${mpi_wrapper_include_text}")

set(openfpm_gpu_example_include_directories
	"${CMAKE_SOURCE_DIR}/openfpm_pdata/src"
	"${CMAKE_SOURCE_DIR}/openfpm_devices/src"
	"${CMAKE_SOURCE_DIR}/openfpm_vcluster/src"
	"${CMAKE_SOURCE_DIR}/openfpm_data/src"
	"${CMAKE_SOURCE_DIR}/openfpm_io/src"
	"${PROJECT_BINARY_DIR}/config"
	${MPI_CXX_INCLUDE_DIRS} ${MPI_C_INCLUDE_DIRS} ${mpi_wrapper_include_dirs}
	${Boost_INCLUDE_DIRS}
	"${PARMETIS_ROOT}/include" "${METIS_ROOT}/include"
	"${HDF5_ROOT}/include" ${LIBHILBERT_INCLUDE_DIRS} ${Vc_INCLUDE_DIR}
	"${VULKAN_INCLUDE_DIR}")

if (ENABLE_NUMERICS)
	list(APPEND openfpm_gpu_example_include_directories ${PETSC_INCLUDES})
endif()

add_executable(gpu_fstep main.cu)
set_source_files_properties(main.cu PROPERTIES LANGUAGE CXX COMPILE_FLAGS "-x c++")
target_compile_features(gpu_fstep PRIVATE cxx_std_17)
target_compile_definitions(gpu_fstep PRIVATE __NVCC__ CUDIFY_USE_METAL CUDA_GPU)
target_include_directories(gpu_fstep BEFORE PRIVATE
	${openfpm_gpu_example_include_directories})
target_link_libraries(gpu_fstep PRIVATE binary_config_ ofpm_pdata vcluster
	ofpmmemory "-framework Foundation" "-framework QuartzCore"
	"-framework Metal" "-framework IOSurface")

include("${CMAKE_SOURCE_DIR}/openfpm_devices/cmake/OpenFPMMoltenVK.cmake")
openfpm_add_moltenvk_kernels(
	TARGET gpu_fstep
	SOURCE main.cu
	HIP_CLANG "${OPENFPM_MOLTENVK_HIP_CLANG}"
	LLVM_OPT "${OPENFPM_MOLTENVK_LLVM_OPT}"
	LLVM_LINK "${OPENFPM_MOLTENVK_LLVM_LINK}"
	CLSPV "${OPENFPM_MOLTENVK_CLSPV}"
	SPIRV_VAL "${OPENFPM_MOLTENVK_SPIRV_VAL}"
	ABI_TOOL "${OPENFPM_MOLTENVK_ABI_TOOL}"
	CHIPSTAR_ROOT "${OPENFPM_CHIPSTAR_ROOT}"
	CHIPSTAR_GENERATED_INCLUDE "${OPENFPM_CHIPSTAR_GENERATED_INCLUDE}"
	COMPILE_DEFINITIONS __NVCC__ CUDIFY_USE_METAL CUDA_GPU
	INCLUDE_DIRECTORIES ${openfpm_gpu_example_include_directories})

set_property(TARGET gpu_fstep PROPERTY BUILD_WITH_INSTALL_RPATH TRUE)
set_property(TARGET gpu_fstep PROPERTY INSTALL_RPATH
	"${METIS_ROOT}/lib;${BOOST_ROOT}/lib;${MPI_ROOT}/lib")
add_custom_command(TARGET gpu_fstep POST_BUILD
	COMMAND install_name_tool -change libmetis.dylib
		@rpath/libmetis.dylib $<TARGET_FILE:gpu_fstep>)

add_test(NAME vector_gpu_first_step_metal
	COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
		$<TARGET_FILE:gpu_fstep>)
