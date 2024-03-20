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
* GPU execution backends include [CUDA](https://developer.nvidia.com/cuda-toolkit), [HIP](https://rocm.docs.amd.com/projects/HIP/en/latest/), [OpenMP](https://www.openmp.org/), [alpaka](https://alpaka.readthedocs.io/en/latest/)
* ... and many others

## Installation
We support MacOS, Linux and Windows subsystem for Linux.

To install please refer to the website instructions on [how to build from source](http://openfpm.mpi-cbg.de/building/)

# Examples and documentation

Examples and documentation are available on the [website](http://openfpm.mpi-cbg.de/news/), in the [online Doxygen documentation](http://ppmcore.mpi-cbg.de/doxygen/openfpm/index.html) and under the folder `example`. 

Example codes include codes for discrete item-based (e.g. Lennard-Jones molecular dynamics) and continuous time and/or space simulations (e.g. Reaction-diffusion-advection equations):

* Dam-break simulation of weakly compressible Navier-Stokes equations in SPH formulation
* Diffusive heat conduction using sparse-grid level set formulation
* 3D Active Fluid simulation
* Gray-Scott reaction-system in 3D
* Hybrid particle-mesh Vortex Method to solve incompressible Navier-Stokes equations
* ... and others

