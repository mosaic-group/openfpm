## OpenFPM: A scalable open-source framework for particle and particle-mesh codes on parallel computers
![OpenFPM.png](OpenFPM.png)

OpenFPM is an open-source software library that facilitates implementing scalable particle and hybrid particle-mesh simulation codes on heterogeneous shared-memory and distributed-memory parallel computer systems.

The library features:

* Scalable serial and parallel data structures for heterogeneous computing systems available CPU and GPU-accelerated hardware
* Operators for particle methods linear differential discretization, e.g. DC-PSE, SPH
* Particle-mesh and mesh-particle interpolation schemes
* Data structures for efficient particle methods simulations, e.g. Cell-List, Verlet-List
* Sparse grid on CPU and GPU
* Support for [PETSc](https://petsc.org/), [Eigen](https://eigen.tuxfamily.org/index.php) linear algebra backends
* Support for ODE integration operators with [Boost.Numeric.Odeint](https://www.boost.org/doc/libs/1_82_0/libs/numeric/odeint/doc/html/index.html)
* Level-set formulation with [Algoim](https://algoim.github.io)
* GPU execution backends include [CUDA](https://developer.nvidia.com/cuda-toolkit), [HIP](https://rocm.docs.amd.com/projects/HIP/en/latest/), [Metal](https://developer.apple.com/metal/) through MoltenVK/SPIR-V, [OpenMP](https://www.openmp.org/), and [alpaka](https://alpaka.readthedocs.io/en/latest/)
* ... and many others

## Installation
We support MacOS, Linux and Windows subsystem for Linux.

To install please refer to the website instructions on [how to build from source](http://openfpm.mpi-cbg.de/building/)

### Selecting the GPU backend

Select the backend when configuring OpenFPM with `CUDA_ON_BACKEND`. The
selection is recorded in the installed CMake package and is used automatically
when downstream projects compile CUDA-style `.cu` sources.

| `CUDA_ON_BACKEND` | Downstream `.cu` compilation |
| --- | --- |
| `CUDA` | NVIDIA CUDA compiler |
| `HIP` | AMD HIP compiler |
| `METAL` | Host code as C++; kernels translated to SPIR-V and executed through MoltenVK |
| `OpenMP` | Cudified source as OpenMP C++ |
| `SEQUENTIAL` | Cudified source as sequential C++ |
| `NONE` | Ordinary C++ with GPU-only guarded code disabled |

Use a separate build and installation directory for each backend. For example,
after providing the dependency paths described in the full installation guide:

```sh
cmake -S . -B build/openfpm-cuda \
  -DCMAKE_INSTALL_PREFIX="$PWD/install/openfpm-cuda" \
  -DCUDA_ON_BACKEND=CUDA \
  -DENABLE_NUMERICS=ON
cmake --build build/openfpm-cuda --parallel
cmake --install build/openfpm-cuda
```

For Metal on macOS, configure the MoltenVK, Vulkan-Headers, clspv/LLVM, and
chipStar paths as part of the OpenFPM build:

```sh
cmake -S . -B build/openfpm-metal \
  -DCMAKE_INSTALL_PREFIX="$PWD/install/openfpm-metal" \
  -DCUDA_ON_BACKEND=METAL \
  -DMOLTENVK_ROOT=/path/to/molten-vk \
  -DVULKAN_HEADERS_ROOT=/path/to/vulkan-headers \
  -DOPENFPM_MOLTENVK_CLSPV=/path/to/clspv \
  -DOPENFPM_MOLTENVK_HIP_CLANG=/path/to/clang++ \
  -DOPENFPM_MOLTENVK_LLVM_OPT=/path/to/opt \
  -DOPENFPM_MOLTENVK_LLVM_LINK=/path/to/llvm-link \
  -DOPENFPM_CHIPSTAR_ROOT=/path/to/chipStar \
  -DENABLE_NUMERICS=ON
cmake --build build/openfpm-metal --parallel
cmake --install build/openfpm-metal
```

## Examples and documentation

Examples and documentation are available on the [website](http://openfpm.mpi-cbg.de/news/), in the [online Doxygen documentation](http://ppmcore.mpi-cbg.de/doxygen/openfpm/index.html) and under the folder `example`. 

### Compiling an example with CMake

Each example directory is a standalone CMake project. Point it at one installed
OpenFPM prefix, then configure and build it out of source:

```sh
cmake -S example/Vector/7_SPH_dlb_gpu -B build/sph-dlb \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PWD/install/openfpm-metal"
cmake --build build/sph-dlb --parallel
cmake --build build/sph-dlb --target run
```

The `run` target launches the executable with the MPI implementation and
runtime library paths recorded by the OpenFPM package, so it does not require
`source openfpm_vars`. It uses one MPI process by default. Select another
process count while configuring the example, for example
`-DOPENFPM_MPI_PROCESSES=4`.

You can use `-Dopenfpm_DIR=/path/to/openfpm/install/cmake` instead of
`CMAKE_PREFIX_PATH`. To switch the example from Metal to CUDA, HIP, or a CPU
backend, point it at an OpenFPM installation configured for that backend; do
not change the example source or add backend-specific targets.

A GPU example uses the same CMake project for every backend:

```cmake
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(my_gpu_example LANGUAGES CXX)

find_package(openfpm CONFIG REQUIRED)

openfpm_add_gpu_executable(my_gpu_example
    SOURCES main.cu)
target_link_libraries(my_gpu_example PRIVATE openfpm::binary_config)
openfpm_add_mpi_run_target(my_gpu_example TARGET_NAME run)
```

Examples that require PETSc or other OpenFPM numerics support request the
component explicitly and link `openfpm::numerics`:

```cmake
find_package(openfpm CONFIG REQUIRED COMPONENTS numerics)
target_link_libraries(my_example PRIVATE openfpm::numerics)
```

CMake reports why the `numerics` component is unavailable when OpenFPM was
built without it or its configured PETSc installation cannot be found. More
example-specific notes are in [`example/README.md`](example/README.md).

Example codes include codes for discrete item-based (e.g. Lennard-Jones molecular dynamics) and continuous time and/or space simulations (e.g. Reaction-diffusion-advection equations):

* Dam-break simulation of weakly compressible Navier-Stokes equations in SPH formulation
* Diffusive heat conduction using sparse-grid level set formulation
* 3D Active Fluid simulation
* Gray-Scott reaction-system in 3D
* Hybrid particle-mesh Vortex Method to solve incompressible Navier-Stokes equations
* ... and others
