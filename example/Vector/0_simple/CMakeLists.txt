cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

cmake_policy(SET CMP0074 NEW)

set(CMAKE_SKIP_BUILD_RPATH TRUE)   ###### <--------- This is absolutely necessary if you use linear algebra

find_package(OpenMP)
find_package(Threads)
find_package(MPI)
find_package(Boost 1.75.0 COMPONENTS program_options iostreams filesystem fiber context)

# CMAKE_PREFIX_PATH can be used to control where to find OpenFPM, the other way is set the path with openfpm_ROOT and activate policy CMP0074
find_package(openfpm 4.1.0 REQUIRED)

if (openfpm_FOUND)
        add_executable(vect  main.cpp)
        target_link_libraries(vect PUBLIC openfpm::binary_config)

        if (OpenMP_FOUND)
                target_link_libraries(vect PUBLIC OpenMP::OpenMP_CXX)
        endif()

        # or

        #target_include_directories(vect PUBLIC ${OPENFPM_INCLUDES})
        #target_compile_definitions(vect PUBLIC ${OPENFPM_DEFINITION})
        #target_link_libraries(vect PUBLIC ${OPENFPM_LIBS})
endif()

