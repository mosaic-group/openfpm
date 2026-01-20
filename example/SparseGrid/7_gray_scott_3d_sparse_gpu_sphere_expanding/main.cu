#define SYNC_BEFORE_TAKE_TIME
#include "Decomposition/Distribution/BoxDistribution.hpp"
#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

/*!
 *
 * \page Grid_3_gs_3D_sparse_gpu_cs_exp Gray Scott in 3D using sparse grids on GPU in expanding complex geometry
 *
 * [TOC]
 *
 * # Solving a gray scott-system in 3D using Sparse grids# {#e3_gs_gray_scott}
 *
 * This example show how to solve a Gray-Scott system in 3D using sparse grids on gpu with complex geometry
 *
 * In figure is the final solution of the problem
 *
 * \htmlonly
<table border="1" bgcolor="black">
  <tr>
    <td>
      <img src="http://ppmcore.mpi-cbg.de/web/images/examples/1_gray_scott_3d_sparse_cs/gs_3d_sparse_cs_section.png" style="width: 500px;" />
    </td>
    <td>
      <img src="http://ppmcore.mpi-cbg.de/web/images/examples/1_gray_scott_3d_sparse_cs/gs_3d_sparse_cs.png" style="width: 500px;" />
    </td>
  </tr>
</table>
\endhtmlonly
 *
 * More or less this example is the same of \ref e3_gs_gray_scott_cs on gpu using what we learned in \ref e3_gs_gray_scott_gpu
 *
 *
 */

#ifdef __NVCC__

constexpr int U = 0;
constexpr int V = 1;
constexpr int U_next = 2;
constexpr int V_next = 3;

typedef CartDecomposition<3,double, CudaMemory, memory_traits_inte, BoxDistribution<3,double> > Dec;

typedef sgrid_dist_id_gpu<3,double,aggregate<double,double,double,double>, CudaMemory,Dec > sgrid_type;

void draw_oscillation_shock(sgrid_type & grid, Box<3,double> & domain)
{
	auto it = grid.getGridIterator();
	Point<3,double> p({1.25,1.25,1.25});

	
//	Point<3,double> u({1.0,0.0,0.0});
//	Box<3,double> channel_box(p3,p1);

	double spacing_x = grid.spacing(0);
	double spacing_y = grid.spacing(1);
	double spacing_z = grid.spacing(2);

	typedef typename GetAddBlockType<sgrid_type>::type InsertBlockT;

	// Draw a shock expanding from 0.4 to 0.8 and than contracting from 0.8 to 0.4
	for (int i = 0 ; i < 100 ; i++)
	{
		Sphere<3,double> sph(p,0.2 + (double)i/160.0);
		Sphere<3,double> sph2(p,0.4 + (double)i/160.0);

		Box<3,size_t> bx;

		for (int j = 0 ; j < 3 ; j++)
		{
			bx.setLow(j,(size_t)((sph.center(j) - 0.4 - (double)i/160.0)/grid.spacing(j)));
			bx.setHigh(j,(size_t)((sph.center(j) + 0.4 + (double)i/160.0)/grid.spacing(j)));
		}

		timer t_add;
		t_add.start();

		grid.addPoints(bx.getKP1(),bx.getKP2(),[spacing_x,spacing_y,spacing_z,sph,sph2] __device__ (int i, int j, int k)
                                {
                                                Point<3,double> pc({i*spacing_x,j*spacing_y,k*spacing_z});

						// Check if the point is in the domain
                                		if (sph2.isInside(pc) )
                                		{
							if (sph.isInside(pc) == false)
							{return true;}
						}

                                                return false;
                                },
                                [] __device__ (InsertBlockT & data, int i, int j, int k)
                                {
                                        data.template get<U>() = 1.0;
                                        data.template get<V>() = 0.0;
                                }
                                );

		t_add.stop();

		timer t_flush;
                t_flush.start();
                grid.template flush<smax_<U>,smax_<V>>(flush_type::FLUSH_ON_DEVICE);
                t_flush.stop();

                timer t_ghost;
                t_ghost.start();
                grid.template ghost_get<U,V>(RUN_ON_DEVICE);
                t_ghost.stop();
                timer t_ghost2;
                t_ghost2.start();
                grid.template ghost_get<U,V>(RUN_ON_DEVICE | SKIP_LABELLING);
                t_ghost2.stop();
                std::cout << t_ghost.getwct() << std::endl;

                std::cout << "TIME ghost1: " << t_ghost.getwct() << "  ghost2: " << t_ghost2.getwct()  << " flush: " <<  t_flush.getwct() << " " << std::endl;


		grid.removeUnusedBuffers();

	}

	std::cout << "Second Pass" <<std::endl;

	for (int i = 0 ; i < 100 ; i++)
	{
		Sphere<3,double> sph(p,0.2 + (double)i/160.0);
		Sphere<3,double> sph2(p,0.4 + (double)i/160.0);

		Box<3,size_t> bx;

		for (int j = 0 ; j < 3 ; j++)
		{
			bx.setLow(j,(size_t)((sph.center(j) - 0.4 - (double)i/160.0)/grid.spacing(j)));
			bx.setHigh(j,(size_t)((sph.center(j) + 0.4 + (double)i/160.0)/grid.spacing(j)));
		}

		timer t_add;
		t_add.start();

		grid.addPoints(bx.getKP1(),bx.getKP2(),[spacing_x,spacing_y,spacing_z,sph,sph2] __device__ (int i, int j, int k)
                                {
                                                Point<3,double> pc({i*spacing_x,j*spacing_y,k*spacing_z});

						// Check if the point is in the domain
                                		if (sph2.isInside(pc) )
                                		{
							if (sph.isInside(pc) == false)
							{return true;}
						}

                                                return false;
                                },
                                [] __device__ (InsertBlockT & data, int i, int j, int k)
                                {
                                        data.template get<U>() = 1.0;
                                        data.template get<V>() = 0.0;
                                }
                                );

		t_add.stop();


		timer t_flush;
                t_flush.start();
		grid.template flush<smax_<U>,smax_<V>>(flush_type::FLUSH_ON_DEVICE);
		t_flush.stop();
//		grid.removeUnusedBuffers();


		timer t_ghost;
		t_ghost.start();
		grid.template ghost_get<U,V>(RUN_ON_DEVICE);
		t_ghost.stop();
		timer t_ghost2;
                t_ghost2.start();
                grid.template ghost_get<U,V>(RUN_ON_DEVICE | SKIP_LABELLING);
                t_ghost2.stop();

		std::cout << "TIME ghost1: " << t_ghost.getwct() << "  ghost2: " << t_ghost2.getwct()  << " flush: " <<  t_flush.getwct() << " " << std::endl;

//		if (i % 10 == 0)
//		{
//			grid.template deviceToHost<U,V>();
//        		grid.write_frame("Final",i);
//		}
	}

}


int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	// domain
	Box<3,double> domain({0.0,0.0,0.0},{2.5,2.5,2.5});
	
	// grid size
        size_t sz[3] = {384,384,384};

	// Define periodicity of the grid
	periodicity<3> bc = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
	
	// Ghost in grid unit
	Ghost<3,long int> g(1);
	
	// deltaT
	double deltaT = 0.025;

	// Diffusion constant for specie U
	double du = 2*1e-5;

	// Diffusion constant for specie V
	double dv = 1*1e-5;

#ifdef TEST_RUN
        // Number of timesteps
        size_t timeSteps = 300;
#else
	// Number of timesteps
        size_t timeSteps = 50000;
#endif

	// K and F (Physical constant in the equation)
        double K = 0.053;
        double F = 0.014;

	grid_sm<3,void> gv({3,1,1});

	sgrid_type grid(sz,domain,g,bc,0,gv);

	grid.template setBackgroundValue<0>(-0.5);
	grid.template setBackgroundValue<1>(-0.5);
	grid.template setBackgroundValue<2>(-0.5);
	grid.template setBackgroundValue<3>(-0.5);
	
	// spacing of the grid on x and y
	double spacing[3] = {grid.spacing(0),grid.spacing(1),grid.spacing(2)};

	draw_oscillation_shock(grid,domain);

	grid.template deviceToHost<U,V>();
	grid.write("Final");

	//! \cond [time stepping] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse_gpu_cs_exp Gray Scott in 3D using sparse grids on GPU in expanding complex geometry
	 *
	 * ## Finalize ##
	 *
	 * Deinitialize the library
	 *
	 * \snippet  SparseGrid/1_gray_scott_3d_sparse_gpu_cs/main.cu finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse_gpu_cs_exp Gray Scott in 3D using sparse grids on GPU in expanding complex geometry
	 *
	 * # Full code # {#code}
	 *
	 * \include SparseGrid/1_gray_scott_3d_sparse_gpu_cs/main.cu
	 *
	 */
}

#else

int main(int argc, char* argv[])
{
        return 0;
}

#endif

