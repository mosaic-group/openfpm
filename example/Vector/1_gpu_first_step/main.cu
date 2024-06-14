/*! \page Vector_1_gpu_first_step Vector 1 GPU first step
 *
 *
 * [TOC]
 *
 *
 * # GPU first steps # {#GPU_first_steps}
 *
 *
 * This example shows how to use GPUs in OpenFPM step by step. To start to use GPUs with vectors_dist, this example is a good
 * starting point. On the other hand we suggest to read the example \ref{simple_vector_example} before this example.
 *
 * ## Data structures {#GPU_data_structures}
 *
 * While a cpu data-structure can be created with vector_dist a gpu data-structure can be created with vector_dist_gpu.
 * The GPU vector_dist expose the same CPU interface with additional functionalities. This means that a vector_dist can be
 * changed to vector_dist_gpu without changing a single line of code. This is an important feature because give us the possibility
 * to change our code from CPU to GPU incrementally step-by-step. A small sections of code can be moved to GPU leaving the rest
 * unchanged. The documentation of vector_dist_gpu is the same ad vector_dist, the extended functionality is documented in
 * vector_dist. Every file containing vector_dist_gpu must be compiled with nvcc compiler, so we suggest to use the extension
 * *.cu for such files and add a rule to compile cu files using nvcc. An example of how to do it can be seen checking the Makefile
 * in this example.
 *
 * While changing from vector_dist to vector_dist_gpu, seem not producing any effect. There are some underling change that take
 * effect:
 *
 * * The internal memory used to allocate memory is not anymore simple heap memory, but is CUDA pinned host memory.
 * * Buffers are allocated for each property.
 * * Vector properties like float[3] and float[3][3] has a different layout in memory from the standard (It will be
 *   clear how and why later in the example)
 *
 * This code snippet below shows how vector_dist_gpu can be used like vector_dist in \ref{simple_vector_example}.
 * In short we create a set of 100 particles (vector_dist_gpu) in 2d from {0.0,0.0} to {1.0,1.0}. Particles are
 * randomly placed in such space. The final map redistribute such particles accordingly to the decomposition.
 *
 * \snippet Vector/1_gpu_first_step/main.cu cpu_like_gpu
 *
 * ### Offload data to GPU ###
 *
 * To offload data to GPU you can use the function hostToDevicePos to offload particle position and hostToDeviceProp to offload
 * properties data. This latest function take template parameters to specify which properties to offload. Here we offload all
 *  the properties (scalar,vector,tensor)
 *
 * \snippet Vector/1_gpu_first_step/main.cu offload_pos_prop
 *
 * Once the data has been offload we can launch a kernel to do some computation. All data-structure gpu ready has a function
 * called toKernel that give the possibility to pass the data-structure to the kernel and be used inside the kernel, like is
 *  used on CPU. Lanuching a kernel on cuda require the subdivision of of a loop in workgroups and threads in the workgroups.
 *  OpenFPM provides the function getDomainIteratorGPU to automatically split the domain particle loop. The members wthr and
 *  thr can be used in the <<<...>>> brachet to launch a CUDA kernel.
 *
 * \snippet Vector/1_gpu_first_step/main.cu launch_domain_it
 *
 * The kernel is the definition of a normal CUDA kernel. We use template parameters for every parameter that is passed with toKernel()
 *
 * \note The definition of the arguments toKernel() as template parameter give us the possibility to use the template engine
 *       to do type deduction and avoid to specify the real-type returned by toKernel()
 *
 * The kernel simply shift the particles by 0.05. Set the scalar properties to the sum of x and y of the "old" particle position,
 * set the vector properties to the old particle position, and set the tensor to several combination of x and y "old" particle
 * position
 *
 * \snippet Vector/1_gpu_first_step/main.cu kernel_translate_fill_prop
 *
 * Once the computation is completed we can ask to reoffload the data from device to host and write the results to file.
 *
 * \note Every file writer requires that the data are offloaded on host memory
 *
 * \snippet Vector/1_gpu_first_step/main.cu device_to_host_write
 *
 * \htmlonly
 * <img src="http://openfpm.mpi-cbg.de/web/images/examples/1_gpu_first_step/output.png"/>
 * \endhtmlonly
 *
 * ## map and ghost_get for multi GPU
 *
 * Until here we saw how to move data from host to device, device to host and how to launch a CUDA kernel on off-loaded data.
 * As previously mentioned vector_dist_gpu has the same CPU interface and so provide the standard function map and ghost_get that work
 * on host pinned memory. Because we want to avoid to move data from GPU to host memory. To avoid it we can use map with the option
 * RUN_DEVICE to redistribute the particles directly on GPU, and ghost_get with RUN_DEVICE to fill ghost particles directly on GPU.
 * In the loop below we see how we can use map on a particle set that is already on GPU. In particular we never offload particles on CPU
 * to do map or ghost_get. We use the kernel translate_fill_prop, to translate the particles and update the properties. The only offload
 * happen every 10 time-step to write on file.
 *
 * \snippet Vector/1_gpu_first_step/main.cu map_and_ghost_get_on_gpu
 *
 * ## RDMA on MPI with CUDA
 *
 * Today MPI implementations are able to do RDMA on GPU memory. This in practice mean that Infiniband card can directly read
 * GPU memory transfer over infiniband and write on the other node directly on GPU, without moving the data to system memory.
 * In practice means that MPI calls can work directly on CUDA device pointers. OpenFPM can exploit this feature if MPI is compiled
 * with CUDA support. To check if MPI is compiled with CUDA support use the function \b is_mpi_rdma_cuda_active() \b
 *
 * \snippet Vector/1_gpu_first_step/main.cu performance_rdma
 *
 * It is good to note that in order to work (return true), some condition must be met.
 *
 * * Because at the moment OpenFPM sense OpenMPI CUDA aware implementation we must define the \b OPENMPI \b macro
 *   \snippet Vector/1_gpu_first_step/main.cu using_openmpi
 *
 * * MPI must be compiled with CUDA support (in general installing OpenFPM with -g should attempt to install OpenMPI with CUDA support)
 * 
 * ## Macro CUDA_LAUNCH
 *
 * When we want to launch a kernel "my_kernel" on CUDA we in general use the Nvidia CUDA syntax
 *
 * my_kernel<<<wthr,thr>>>(arguments ... )
 *
 * Where wthr is the number of workgroups and thr is the number of threads in a workgroup and arguments... are the arguments to pass to the kernel.
 * Equivalently we can launch a kernel with the macro CUDA_LAUNCH_DIM3(my_kernel,wthr,thr,arguments...) or CUDA_LAUNCH(my_kernel,ite,arguments) where
 * ite has been taken using getDomainIteratorGPU. There are several advantage on using CUDA_LAUNCH. The first advantage in using the macro is enabling SE_CLASS1
 * all kernel launch become synchronous and an error check is performed before continue to the next kernel making debugging easier. Another feature is the possibility
 * to run CUDA code on CPU without a GPU. compiling with "CUDA_ON_CPU=1 make" (Note openfpm must be compiled with GPU support (-g)  or with CUDA_ON_CPU support
 * (-c "... --enable_cuda_on_cpu"). You can compile this example on CPU. You do not have to change a single line of code for this example. (Check the video to see this
 * feature in action). All the openfpm GPU example and CUDA example can run on CPU if they use CUDA_LAUNCH as macro. We are planning to support
 * AMD GPUs as well using this system.
 *
 * ## Full code ## {#code_e0_sim}
 *
 * \include Vector/1_gpu_first_step/main.cu
 *
 */

#ifdef __NVCC__

//! \cond [using_openmpi] \endcond
#define OPENMPI
//! \cond [using_openmpi] \endcond

//#define EXTERNAL_SET_GPU <----- In case you want to distribute the GPUs differently from the default

#include "Vector/vector_dist.hpp"

//! \cond [kernel_translate_fill_prop] \endcond

template<typename vector_type>
__global__ void translate_fill_prop(vector_type vd)
{
	auto p = GET_PARTICLE(vd);

	vd.template getProp<0>(p) = vd.getPos(p)[0] + vd.getPos(p)[1];

	vd.template getProp<1>(p)[0] = vd.getPos(p)[0];
	vd.template getProp<1>(p)[1] = vd.getPos(p)[1];

	vd.template getProp<2>(p)[0][0] = vd.getPos(p)[0];
	vd.template getProp<2>(p)[0][1] = vd.getPos(p)[1];
	vd.template getProp<2>(p)[1][0] = vd.getPos(p)[0] + vd.getPos(p)[1];
	vd.template getProp<2>(p)[1][1] = vd.getPos(p)[1] - vd.getPos(p)[0];

	vd.getPos(p)[0] += 0.01f;
	vd.getPos(p)[1] += 0.01f;
}

//! \cond [kernel_translate_fill_prop] \endcond

int main(int argc, char* argv[])
{
    // OpenFPM GPU distribution

    // OpenFPM by default select GPU 0 for process 0, gpu 1 for process 1 and so on ... . In case of multi-node is the same each node has
    // has a group of processes and these group of processes are distributed across the available GPU on that node.

    // If you want to override this behaviour use #define EXTERNAL_SET_GPU at the very beginning of the program and use
    // cudaSetDevice to select the GPU for that particular process before openfpm_init
    // Note: To get the process number do MPI_Init and and use the MPI_Comm_rank. VCluster is not available before openfpm_init
    // A code snippet in case we want to skip GPU 0
    // MPI_Init(&argc,&argv);
    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    // cudaSetDevice(1+rank);

	//! \cond [cpu_like_gpu] \endcond

    // initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<2,float> g(0.05);

    vector_dist_gpu<2,float, aggregate<float,float[2],float[2][2]> > vd(100,domain,bc,g);

	// the scalar is the element at position 0 in the aggregate
	const int scalar = 0;

	// the vector is the element at position 1 in the aggregate
	const int vector = 1;

	// the tensor is the element at position 2 in the aggregate
	const int tensor = 2;

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		// we define x, assign a random position between 0.0 and 1.0
		vd.getPos(key)[0] = (float)rand() / RAND_MAX;

		// we define y, assign a random position between 0.0 and 1.0
		vd.getPos(key)[1] = (float)rand() / RAND_MAX;

		// next particle
		++it;
	}

	vd.map();

	//! \cond [cpu_like_gpu] \endcond

	//! \cond [offload_pos_prop] \endcond

	vd.hostToDevicePos();
	vd.template hostToDeviceProp<scalar,vector,tensor>();

	//! \cond [offload_pos_prop] \endcond

	//! \cond [launch_domain_it] \endcond

	auto ite = vd.getDomainIteratorGPU();
	// translate_fill_prop<<<ite.wthr,ite.thr>>>(vd.toKernel());
	CUDA_LAUNCH(translate_fill_prop,ite,vd.toKernel());

	//! \cond [launch_domain_it] \endcond

	//! \cond [device_to_host_write] \endcond

	vd.deviceToHostPos();
	vd.deviceToHostProp<0,1,2>();

	// We write on a file
	vd.write("output");

	//! \cond [device_to_host_write] \endcond

	//! \cond [map_and_ghost_get_on_gpu] \endcond

	for (int j = 0 ; j < 100 ; j++)
	{
		auto ite = vd.getDomainIteratorGPU();
		// translate_fill_prop<<<ite.wthr,ite.thr>>>(vd.toKernel());
		CUDA_LAUNCH(translate_fill_prop,ite,vd.toKernel());

		vd.map(RUN_ON_DEVICE);
		vd.template ghost_get<0,1,2>(RUN_ON_DEVICE);

		if ( j % 10 == 0)
		{
			// offload to host
			vd.deviceToHostPos();
			vd.template deviceToHostProp<0,1,2>();

			// write
			vd.write_frame("output_f",j);
		}
	}

	//! \cond [map_and_ghost_get_on_gpu] \endcond

	//! \cond [performance_rdma] \endcond

	bool active = is_mpi_rdma_cuda_active();

	std::cout << "Is MPI rdma active on CUDA " << active << std::endl;

	//! \cond [performance_rdma] \endcond

	openfpm_finalize();
}

#else

int main(int argc, char* argv[])
{
        return 0;
}

#endif
