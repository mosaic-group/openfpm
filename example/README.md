# Building OpenFPM examples

Each source example is a standalone CMake consumer of the installed OpenFPM
package. Configure and build an example by pointing CMake at the installation:

```sh
cmake -S example/Vector/7_SPH_dlb_gpu -B build/sph-dlb \
  -Dopenfpm_DIR=/path/to/openfpm/install/cmake
cmake --build build/sph-dlb
```

GPU examples always keep their CUDA-style `.cu` source. The installed OpenFPM
configuration selects whether CMake compiles that source with CUDA, HIP,
Metal/SPIR-V, OpenMP, or the sequential backend.

Examples that use the optional numerics module request it explicitly:

```cmake
find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)
target_link_libraries(my_example PRIVATE openfpm::numerics)
```

CMake reports the `numerics` component as unavailable when the installed
OpenFPM was built without numerics or when its configured PETSc installation
is missing. Core-only consumers continue to use `openfpm::binary_config` and
do not require the numerics component.

The two historical libquadmath examples are disabled by default. Configure
their directories with `-DOPENFPM_BUILD_QUADMATH_EXAMPLE=ON` on systems that
provide libquadmath.
