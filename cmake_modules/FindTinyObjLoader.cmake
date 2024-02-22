# - Try to find LibHilbert
# Once done this will define
#
#  TINYOBJLOADER_FOUND        - system has TinyObjLoader
#  TINYOBJLOADER_INCLUDE_DIRS - include directories for TinyObjLoader
#  TINYOBJLOADER_LIBRARY_DIRS - library directories for TinyObjLoader
#  TINYOBJLOADER_LIBRARIES    - libraries for TinyObjLoader
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


if (TINYOBJLOADER_FOUND OR TARGET tinyobjloader)
	return()
endif()

add_library(tinyobjloader INTERFACE IMPORTED)

# Add libraries (static)
set(_libs "-L${TINYOBJLOADER_ROOT}/lib64 -ltinyobjloader")
set_property(TARGET tinyobjloader PROPERTY INTERFACE_LINK_LIBRARIES "${_libs}")


# Create LibHilbert test program
set(TINYOBJLOADER_TEST_LIB_CPP
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/tinyobjloader_test_lib.cpp")


file(WRITE ${TINYOBJLOADER_TEST_LIB_CPP} "
#define TINYOBJLOADER_IMPLEMENTATION
#include \"tiny_obj_loader.h\"


int main()
{
  std::string inputfile = \"null.obj\";
  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;
  
  std::string err;
  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, inputfile.c_str());

  return 0;
}
")

# Try to run test program (static linking)
try_run(
	TINYOBJLOADER_TEST_LIB_EXITCODE
	TINYOBJLOADER_TEST_LIB_COMPILED
      	${CMAKE_CURRENT_BINARY_DIR}
	${TINYOBJLOADER_TEST_LIB_CPP}
      CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES:STRING=${TINYOBJLOADER_ROOT}/include"
      "-DLINK_LIBRARIES:STRING=${TINYOBJLOADER_ROOT}/lib"
      LINK_LIBRARIES tinyobjloader
      COMPILE_OUTPUT_VARIABLE TINYOBJLOADER_TEST_LIB_COMPILE_OUTPUT
      RUN_OUTPUT_VARIABLE TINYOBJLOADER_TEST_LIB_OUTPUT)


if (TINYOBJLOADER_TEST_LIB_COMPILED AND TINYOBJLOADER_TEST_LIB_EXITCODE EQUAL 0)
	message(STATUS "Test TinyObjLoader_TEST_RUNS static linking - Success")
	    set(TINYOBJLOADER_TEST_RUNS TRUE)
	    set(TINYOBJLOADER_FOUND TRUE)
	    set(TINYOBJLOADER_INCLUDE_DIRS ${TINYOBJLOADER_ROOT}/include)
	    set(TINYOBJLOADER_LIBRARY_DIRS ${TINYOBJLOADER_ROOT}/lib64)
	    set(TINYOBJLOADER_LIBRARIES -ltinyobjloader)
else()
	message(STATUS "Test TinyObjLoader_TEST_RUNS static linking - Failed")
	    set(TINYOBJLOADER_TEST_RUNS FALSE)
endif()

