/*! \page Vector_9_gpu_cuda_interop Vector 9 GPU cuda interoperability
 *
 *
 * [TOC]
 *
 *
 * # GPU CUDA inter-operability with plain arrays # {#GPU_9_cuda_interop}
 *
 * OpenFPM provide the possibility to operate with CUDA using plain arrays. In particular we can ask to the distributed
 * data-structure to return a CUDA device pointer to the data. Before operate with such pointer we must understand
 * how vector_dist_gpu store data internally in order to correctly read data from such pointer
 *
 * ## Array striding {#e9_array_stride}
 *
 * To understand how vector_dist_gpu store data, we will print the address in memory of each element. Let start for printing
 *  the address of the first particle for all the properties
 *
 * \snippet Vector/9_gpu_cuda_interop/main.cu first_particle_prop_zero_and_one_two
 *
 * Running the program on one process we get
 *
 * \code
First particle property 0, address: 0x7f8d63c00400
First particle property 1, address: 0x7f8d63c00600
First particle property 2, address: 0x7f8d83400000
 * \endcode
 *
 * As we can see the scalar property, vector and tensor properties are not nearly contiguous, the reason is that every properties use
 * its own CUDA buffer, and each property can be off-loaded separately.
 *
 * Now we check how the component of the vector are stored in memory, printing the address of the components for the vector and for
 * the tensor property.
 *
 * \snippet Vector/9_gpu_cuda_interop/main.cu first_particle_vector_tensor_layout
 *
 * The output that we can obtain is something like
 *
  \code
Capacity internal vector: 128
First particle property 1 component 0, address: 0x7f8d63c00600
First particle property 1 component 1, address: 0x7f8d63c00800
First particle property 2 component 00, address: 0x7f8d83400000
First particle property 2 component 01, address: 0x7f8d83400200
First particle property 2 component 10, address: 0x7f8d83400400
First particle property 2 component 11, address: 0x7f8d83400600
Second particle property 1 component 0, address: 0x7f8d63c00604
Second particle property 1 component 1, address: 0x7f8d63c00804
Second particle property 2 component 00, address: 0x7f8d83400004
Second particle property 2 component 01, address: 0x7f8d83400204
Second particle property 2 component 10, address: 0x7f8d83400404
Second particle property 2 component 11, address: 0x7f8d83400604
  \endcode
 *
 * As we can see the vector property of first particle component y is not contiguous to x, but is 0x200 = 4 byte * 128 offset from
 * the component x. What is contiguous to particle 0 component x is particle 1 component x
 *
 * \note This is in general hidden and transparent to the user. Infact in the example we have shown, we were able to create a distributed vector and
 *       compute on it without know how vector_dist store data. It only become necessary if you want to use CUDA with plain primitive arrays
 *
 * There is a reason why vector_dist_gpu use this layout and is because of memory coalesced access. Suppose you want to access
 * a vector property in the GPU kernel like this
 *
 * \code
 *
 * vd.template getProp<vector>(p)[0]
 *
 * \endcode
 *
 * In general what we do is to map the particle index p to a GPU thread that handle that particle. Doing so let see what happen
 * when one SM hit that instruction using the standard layout.
 *
 * \verbatim
                          Memory                                                              Memory

   particle 0  [0]x        0x000      <------- Access thread 0         particle 0  [0]x        0x000      <------- Access thread 0
               [1]y        0x004                                       particle 1  [0]x        0x004      <------- Access thread 1
   particle 1  [0]x        0x008      <------- Access thread 1         particle 2  [0]x        0x008      <------- Access thread 2
               [1]y          .                                         particle 3  [0]x        0x00C      <------- Access thread 3
   particle 2  [0]x          .                                                  .
               [1]y          .                                                  .
               .             .                                                  .
               .             .                                                  .
               .             .                                                  .
   particle N  [0]x          .       <-------- Access thread N                  .
               [1]y          .                                         particle 0  [1]y

                  Case A                                                            Case B
 * \endverbatim
 *
 * As we can see from the image in case A there is a jump of 4 byte compared of Case B. And this mean that the instruction will read double
 *  of the cache lines compared to case B.
 *
 *  Remain to understand why having 100 particles the component y stay at 4 * 128 = 512 byte instead of 4 * 100 = 400 byte. One power 2 alignment
 *  the other is instead related to the internal preallocated vector buffer. Suppose to have a vector with 4 particles and we want to add one
 *  particle at the end. Because we do not have space in theory we have to create a vector of 5 elements, copy the 4 elements in the new vector
 *  and add the last elements. This is clearly expensive when the vector become big, copy the full vector to just one element would not make sense.
 *  OpenFPM use by default a policy to expand the vector by a factor (default = 2) to guarantee that if a vector with N elements starting from
 *  an a vector of size 0 have cost O(N).
 *
 * \note OpenFPM by default does not operate any attempt to expand the virtual address space of the structure to avoid copy
 *
 * ## Interoperability with CUDA {#e9_interop_cuda}
 *
 * Now that we understood the structure of the device pointer, we can see how we can use the internal device pointer in a cuda kernel.
 * We now launch a kernel just to print the information inside the buffer. To get the device CUDA pointer we can use the combo functions
 * \b getPropVector() \b to the the internal propetries vector follow  by \b getDeviceBuffer<0>() \b that return the CUDA device
 * buffer for the property 0
 *
 * \snippet Vector/9_gpu_cuda_interop/main.cu print_50
 *
 * the kernel print the information of particle 50. To note how we pass primitive arrays to the kernel and we use capacity to
 * access the component of vector and the tensor accordingly  to what we explained in array striding
 *
 * \snippet Vector/9_gpu_cuda_interop/main.cu print_data_kernel
 *
 *
 *
 * ## Full code ## {#code_e9_sim}
 *
 * \include Vector/9_gpu_cuda_interop/main.cu
 *
 */

#ifdef __NVCC__

#include "Vector/vector_dist.hpp"
#include "util/cuda_util.hpp"

//! [print_data_kernel]

__global__ void print_data_particle_50(float * scalar, float * vector, float * tensor, size_t capacity)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p == 50)
	{
		printf("Scalar particle %d = %f\n",p,scalar[p]);

		printf("Vector particle %d = %f\n",p,vector[p]);
		printf("Vector particle %d = %f\n",p,vector[p + capacity]);

		printf("Tensor particle %d = %f\n",p,tensor[p + (0*2 + 0)*capacity]);
		printf("Tensor particle %d = %f\n",p,tensor[p + (0*2 + 1)*capacity]);
		printf("Tensor particle %d = %f\n",p,tensor[p + (1*2 + 0)*capacity]);
		printf("Tensor particle %d = %f\n",p,tensor[p + (1*2 + 1)*capacity]);
	}
}

//! [print_data_kernel]

int main(int argc, char* argv[])
{
    // initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<2,float> g(0.05);

    vector_dist_gpu<2,float, aggregate<float,float[2],float[2][2]> > vd(100,domain,bc,g);

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		// we define x, assign a random position between 0.0 and 1.0
		vd.getPos(p)[0] = (float)rand() / RAND_MAX;

		// we define y, assign a random position between 0.0 and 1.0
		vd.getPos(p)[1] = (float)rand() / RAND_MAX;

		vd.template getProp<0>(p) = vd.getPos(p)[0] + vd.getPos(p)[1];

		vd.template getProp<1>(p)[0] = vd.getPos(p)[0];
		vd.template getProp<1>(p)[1] = vd.getPos(p)[1];

		vd.template getProp<2>(p)[0][0] = vd.getPos(p)[0];
		vd.template getProp<2>(p)[0][1] = vd.getPos(p)[1];
		vd.template getProp<2>(p)[1][0] = vd.getPos(p)[0] + vd.getPos(p)[1];
		vd.template getProp<2>(p)[1][1] = vd.getPos(p)[1] - vd.getPos(p)[0];

		// next particle
		++it;
	}

	vd.map();

	//! \cond [map_and_ghost_get_on_gpu] \endcond

	//! \cond [first_particle_prop_zero_and_one_two] \endcond

	std::cout << "First particle property 0, address: " << &vd.template getProp<0>(0) << std::endl;
	std::cout << "First particle property 1, address: " << &vd.template getProp<1>(0)[0] << std::endl;
	std::cout << "First particle property 2, address: " << &vd.template getProp<2>(0)[0][0] << std::endl;

	//! \cond [first_particle_prop_zero_and_one_two] \endcond

	//! \cond [first_particle_vector_tensor_layout] \endcond

	std::cout << "Capacity internal vector: " << vd.getPropVector().capacity() << std::endl;

	std::cout << "First particle property 1 component 0, address: " << &vd.template getProp<1>(0)[0] << std::endl;
	std::cout << "First particle property 1 component 1, address: " << &vd.template getProp<1>(0)[1] << std::endl;

	std::cout << "First particle property 2 component 00, address: " << &vd.template getProp<2>(0)[0][0] << std::endl;
	std::cout << "First particle property 2 component 01, address: " << &vd.template getProp<2>(0)[0][1] << std::endl;
	std::cout << "First particle property 2 component 10, address: " << &vd.template getProp<2>(0)[1][0] << std::endl;
	std::cout << "First particle property 2 component 11, address: " << &vd.template getProp<2>(0)[1][1] << std::endl;

	std::cout << "Second particle property 1 component 0, address: " << &vd.template getProp<1>(1)[0] << std::endl;
	std::cout << "Second particle property 1 component 1, address: " << &vd.template getProp<1>(1)[1] << std::endl;

	std::cout << "Second particle property 2 component 00, address: " << &vd.template getProp<2>(1)[0][0] << std::endl;
	std::cout << "Second particle property 2 component 01, address: " << &vd.template getProp<2>(1)[0][1] << std::endl;
	std::cout << "Second particle property 2 component 10, address: " << &vd.template getProp<2>(1)[1][0] << std::endl;
	std::cout << "Second particle property 2 component 11, address: " << &vd.template getProp<2>(1)[1][1] << std::endl;

	//! \cond [first_particle_vector_tensor_layout] \endcond

	std::cout << std::endl;

	//! \cond [print_50] \endcond

	vd.template hostToDeviceProp<0,1,2>();

	CUDA_LAUNCH_DIM3(print_data_particle_50,100,1,(float *)vd.getPropVector().template getDeviceBuffer<0>(),
			               (float *)vd.getPropVector().template getDeviceBuffer<1>(),
			               (float *)vd.getPropVector().template getDeviceBuffer<2>(),
			               vd.getPropVector().capacity());

	//! \cond [print_50] \endcond

	openfpm_finalize();
}

#else

int main(int argc, char* argv[])
{
        return 0;
}

#endif
