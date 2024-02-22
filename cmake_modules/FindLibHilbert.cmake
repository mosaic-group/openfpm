# - Try to find LibHilbert
# Once done this will define
#
#  LIBHILBERT_FOUND        - system has LibHilbert
#  LIBHILBERT_INCLUDE_DIRS - include directories for PETSc
#  LIBHILBERT_LIBRARY_DIRS - library directories for PETSc
#  LIBHILBERT_LIBRARIES    - libraries for PETSc
#
#=============================================================================
# Copyright (C) 2010-2016 Pietro Incardona
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#=============================================================================


if (LIBHILBERT_FOUND)
	return()
endif()

add_library(libhilbert INTERFACE IMPORTED)

# Add libraries (static)
set(_libs "-L${LIBHILBERT_ROOT}/lib -llibhilbert")
set_property(TARGET libhilbert PROPERTY INTERFACE_LINK_LIBRARIES "${_libs}")


# Create LibHilbert test program
set(LIBHILBERT_TEST_LIB_CPP
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/libhilbert_test_lib.cpp")


file(WRITE ${LIBHILBERT_TEST_LIB_CPP} "
extern \"C\"
{
#include \"hilbertKey.h\"
}
  int main()
{
  //An integer to handle errors
  int err;

  //Array to handle output
  uint64_t nextCoord[2];

  //Get the coordinates of the next cell
  getIntCoordFromHKey(nextCoord, 4, 2, 0, &err);

  return 0;
}
")

# Try to run test program (static linking)
try_run(
	LIBHILBERT_TEST_LIB_EXITCODE
	LIBHILBERT_TEST_LIB_COMPILED
      	${CMAKE_CURRENT_BINARY_DIR}
        ${LIBHILBERT_TEST_LIB_CPP}
      CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES:STRING=${LIBHILBERT_ROOT}/include"
      "-DLINK_LIBRARIES:STRING=${LIBHILBERT_ROOT}/lib"
      LINK_LIBRARIES libhilbert
      COMPILE_OUTPUT_VARIABLE LIBHILBERT_TEST_LIB_COMPILE_OUTPUT
      RUN_OUTPUT_VARIABLE LIBHILBERT_TEST_LIB_OUTPUT)

if (NOT LIBHILBERT_TEST_LIB_COMPILED)
	message("Error in compiling libhilbert test: ${LIBHILBERT_TEST_LIB_COMPILE_OUTPUT}")
endif()

if (LIBHILBERT_TEST_LIB_COMPILED AND LIBHILBERT_TEST_LIB_EXITCODE EQUAL 0)
	    message(STATUS "Test LibHilbert_TEST_RUNS static linking - Success")
	    set(LIBHILBERT_TEST_RUNS TRUE)
	    set(LIBHILBERT_FOUND TRUE)
	    set(LIBHILBERT_INCLUDE_DIRS ${LIBHILBERT_ROOT}/include)
	    set(LIBHILBERT_LIBRARY_DIRS ${LIBHILBERT_ROOT}/lib)
	    set(LIBHILBERT_LIBRARIES -llibhilbert)
else()
	    message(STATUS "Test LibHilbert_TEST_RUNS static linking - Failed")
	    set(LIBHILBERT_TEST_RUNS FALSE)
endif()

