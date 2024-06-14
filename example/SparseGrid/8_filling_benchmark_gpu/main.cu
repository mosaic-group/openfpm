#define VCLUSTER_PERF_REPORT
#define SYNC_BEFORE_TAKE_TIME
#define ENABLE_GRID_DIST_ID_PERF_STATS
#include "Decomposition/Distribution/BoxDistribution.hpp"
#include "util/cuda_util.hpp"
#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

/*!
 *
 * \page Grid_3_gs_3D_sparse_gpu Gray Scott in 3D using sparse grids on GPU
 *
 * [TOC]
 *
 * # Solving a gray scott-system in 3D using Sparse grids on gpu # {#e3_gs_gray_scott_gpu}
 *
 * This example show how to solve a Gray-Scott system in 3D using sparse grids on gpu
 *
 * In figure is the final solution of the problem
 *
 * \htmlonly
 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/gray_scott_3d/gs_alpha.png"/>
 * \endhtmlonly
 *
 * More or less this example is the adaptation of the dense example in 3D
 *
 * \see \ref Grid_3_gs_3D
 *
 * # Initializetion
 *
 * On gpu we can add points using the function addPoints this function take 2 lamda functions the first take 3 arguments (in 3D)
 * i,j,k these are the global coordinates for a point. We can return either true either false. In case of true the point is
 * created in case of false the point is not inserted. The second lamda is instead used to initialize the point inserted.
 * The arguments of the second lambda are the data argument we use to initialize the point and the global coordinates i,j,k
 *
 * After we add the points we have to flush the added points. This us achieved using the function flush the template parameters
 * indicate how we have to act on the points. Consider infact we are adding points already exist ... do we have to add it using the max
 * or the min. **FLUSH_ON_DEVICE** say instead that the operation is performed using the GPU
 *
 * \snippet SparseGrid/1_gray_scott_3d_sparse_gpu/main.cu create points
 *
 * The function can also called with a specified range
 *
 * \snippet SparseGrid/1_gray_scott_3d_sparse_gpu/main.cu create points sub
 *
 * # Update
 *
 * to calculate the right-hand-side we use the function **conv2** this function can be used to do a convolution that involve
 * two properties
 *
 * The function accept a lambda function where the first 2 arguments are the output of the same type of the two property choosen.
 *
 * The arguments 3 and 4 contain the properties of two selected properties. while i,j,k are the coordinates we have to calculate the
 * convolution. The call **conv2** also accept template parameters the first two indicate the source porperties, the other two are the destination properties. While the
 * last is the extension of the stencil. In this case we use 1.
 *
 * The lambda function is defined as
 *
 * \snippet SparseGrid/1_gray_scott_3d_sparse_gpu/main.cu lambda
 *
 * and used in the body loop
 *
 * \snippet SparseGrid/1_gray_scott_3d_sparse_gpu/main.cu body
 *
 */

#ifdef __NVCC__

constexpr int U = 0;
constexpr int V = 1;

constexpr int U_next = 2;
constexpr int V_next = 3;

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

typedef CartDecomposition<3,float, CudaMemory, memory_traits_inte, BoxDistribution<3,float> > Dec;

typedef sgrid_dist_id_gpu<3,float,aggregate<float>,CudaMemory, Dec> SparseGridType;

void init(SparseGridType & grid, Box<3,float> & domain)
{
	//! \cond [create points] \endcond

	typedef typename GetAddBlockType<SparseGridType>::type InsertBlockT;

	for (int i = 0 ; i < 10 ; i++)
	{
		timer t;
		t.start();

		grid.addPoints([] __device__ (int i, int j, int k)
			        {
						return true;
			        },
			        [] __device__ (InsertBlockT & data, int i, int j, int k)
			        {
			        	data.template get<U>() = 1.0;
			        }
			        );


		grid.template flush<smax_<U>>(flush_type::FLUSH_ON_DEVICE);

		t.stop();

		std::cout << "Time populate: " << t.getwct()  << std::endl;

        	timer t2;
		cudaDeviceSynchronize();
        	t2.start();

        	grid.addPoints([] __device__ (int i, int j, int k)
                                {
                                                return true;
                                },
                                [] __device__ (InsertBlockT & data, int i, int j, int k)
                                {
                                        data.template get<U>() = 5.0;
                                }
                                );


        	grid.template flush<sRight_<U>>(flush_type::FLUSH_ON_DEVICE);

		t2.stop();

		std::cout << "Time populate: " << t2.getwct()  << std::endl;
	}
}


int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	// domain
	Box<3,float> domain({0.0,0.0,0.0},{2.5,2.5,2.5});
	
	// grid size
        size_t sz[3] = {512,512,512};

	// Define periodicity of the grid
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	
	// Ghost in grid unit
	Ghost<3,long int> g(1);
	
	// deltaT
	float deltaT = 0.25;

	// Diffusion constant for specie U
	float du = 2*1e-5;

	// Diffusion constant for specie V
	float dv = 1*1e-5;

	// Number of timesteps
#ifdef TEST_RUN
	size_t timeSteps = 300;
#else
        size_t timeSteps = 15000;
#endif

	// K and F (Physical constant in the equation)
    float K = 0.053;
    float F = 0.014;

	SparseGridType grid(sz,domain,g,bc);

	// spacing of the grid on x and y
	float spacing[3] = {grid.spacing(0),grid.spacing(1),grid.spacing(2)};

	init(grid,domain);

	// sync the ghost
	grid.deviceToHost<U>();
	grid.write("final");
	grid.print_stats();

	//! \cond [time stepping] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse Gray Scott in 3D
	 *
	 * ## Finalize ##
	 *
	 * Deinitialize the library
	 *
	 * \snippet Grid/3_gray_scott/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse Gray Scott in 3D
	 *
	 * # Full code # {#code}
	 *
	 * \include Grid/3_gray_scott_3d/main.cpp
	 *
	 */
}

#else

int main(int argc, char* argv[])
{
        return 0;
}

#endif

