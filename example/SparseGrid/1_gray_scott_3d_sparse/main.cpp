 /*! \page SparseGrid SparseGrid
 *
 * \subpage Grid_3_gs_3D_sparse
 * \subpage Grid_3_gs_3D_sparse_gpu
 * \subpage Grid_3_gs_3D_sparse_opt
 * \subpage Grid_3_gs_3D_sparse_gpu_opt
 * \subpage Grid_3_gs_3D_sparse_cs
 * \subpage Grid_3_gs_3D_sparse_gpu_cs
 * \subpage Grid_3_gs_3D_sparse_gpu_cs_opt
 * \subpage Grid_3_gs_3D_sparse_gpu_cs_exp
 * \subpage inhomogeneous_diffusion_porous_catalyst_CaCO3
 * \subpage heat_conduction_RPC
 *
 */

#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

/*!
 *
 * \page Grid_3_gs_3D_sparse Gray Scott in 3D using sparse grids on CPU
 *
 * [TOC]
 *
 * # Solving a gray scott-system in 3D using Sparse grids on GPU optimized # {#1_gray_scott_3d_sparse:}
 *
 * This example show how to solve a Gray-Scott system in 3D using sparse grids on gpu (optimized). The problem is the same as \ref Grid_3_gs_3D
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
 *
 * 
 *
 *
 * 
 *
 */

constexpr int U = 0;
constexpr int V = 1;

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

void init(sgrid_dist_id<3,double,aggregate<double,double> > & Old, sgrid_dist_id<3,double,aggregate<double,double> > & New, Box<3,double> & domain)
{
	auto it = Old.getGridIterator();

	while (it.isNext())
	{
		// Get the local grid key
		auto key = it.get_dist();

		// Old values U and V
		Old.template insert<U>(key) = 1.0;
		Old.template insert<V>(key) = 0.0;

		// Old values U and V
		New.template insert<U>(key) = 0.0;
		New.template insert<V>(key) = 0.0;

		++it;
	}

	long int x_start = Old.size(0)*1.55f/domain.getHigh(0);
	long int y_start = Old.size(1)*1.55f/domain.getHigh(1);
	long int z_start = Old.size(1)*1.55f/domain.getHigh(2);

	long int x_stop = Old.size(0)*1.85f/domain.getHigh(0);
	long int y_stop = Old.size(1)*1.85f/domain.getHigh(1);
	long int z_stop = Old.size(1)*1.85f/domain.getHigh(2);

	grid_key_dx<3> start({x_start,y_start,z_start});
	grid_key_dx<3> stop ({x_stop,y_stop,z_stop});
	auto it_init = Old.getGridIterator(start,stop);

	while (it_init.isNext())
	{
		auto key = it_init.get_dist();

                Old.template insert<U>(key) = 0.5 + (((double)std::rand())/RAND_MAX -0.5)/10.0;
                Old.template insert<V>(key) = 0.25 + (((double)std::rand())/RAND_MAX -0.5)/20.0;

		++it_init;
	}
}


int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	// domain
	Box<3,double> domain({0.0,0.0,0.0},{2.5,2.5,2.5});
	
	// grid size
        size_t sz[3] = {128,128,128};

	// Define periodicity of the grid
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	
	// Ghost in grid unit
	Ghost<3,long int> g(1);
	
	// deltaT
	double deltaT = 1;

	// Diffusion constant for specie U
	double du = 2*1e-5;

	// Diffusion constant for specie V
	double dv = 1*1e-5;

	// Number of timesteps
#ifdef TEST_RUN
	size_t timeSteps = 200;
#else
        size_t timeSteps = 5000;
#endif

	// K and F (Physical constant in the equation)
        double K = 0.053;
        double F = 0.014;

	sgrid_dist_id<3, double, aggregate<double,double>> Old(sz,domain,g,bc);

	// New grid with the decomposition of the old grid
        sgrid_dist_id<3, double, aggregate<double,double>> New(Old.getDecomposition(),sz,g);

	
	// spacing of the grid on x and y
	double spacing[3] = {Old.spacing(0),Old.spacing(1),Old.spacing(2)};

	init(Old,New,domain);

	// sync the ghost
	size_t count = 0;
	Old.template ghost_get<U,V>();

	// because we assume that spacing[x] == spacing[y] we use formula 2
	// and we calculate the prefactor of Eq 2
	double uFactor = deltaT * du/(spacing[x]*spacing[x]);
	double vFactor = deltaT * dv/(spacing[x]*spacing[x]);

	Old.write("Init_condition");

	timer tot_sim;
	tot_sim.start();

	auto & v_cl = create_vcluster(); 

	for (size_t i = 0; i < timeSteps; ++i)
	{
		if (v_cl.rank() == 0)
		{std::cout << "STEP: " << i << std::endl;}
/*		if (i % 300 == 0)
		{
			std::cout << "STEP: " << i << std::endl;
			Old.write_frame("out",i);
		}*/

		//! \cond [stencil get and use] \endcond

		auto it = Old.getDomainIterator();

		while (it.isNext())
		{
			// center point
			auto Cp = it.get();

			// plus,minus X,Y,Z
			auto mx = Cp.move(0,-1);
			auto px = Cp.move(0,+1);
			auto my = Cp.move(1,-1);
			auto py = Cp.move(1,1);
			auto mz = Cp.move(2,-1);
			auto pz = Cp.move(2,1);

			// update based on Eq 2
			New.insert<U>(Cp) = Old.get<U>(Cp) + uFactor * (
										Old.get<U>(mz) +
										Old.get<U>(pz) +
										Old.get<U>(my) +
										Old.get<U>(py) +
										Old.get<U>(mx) +
										Old.get<U>(px) -
										6.0*Old.get<U>(Cp)) +
										- deltaT * Old.get<U>(Cp) * Old.get<V>(Cp) * Old.get<V>(Cp) +
										- deltaT * F * (Old.get<U>(Cp) - 1.0);


			// update based on Eq 2
			New.insert<V>(Cp) = Old.get<V>(Cp) + vFactor * (
										Old.get<V>(mz) +
										Old.get<V>(pz) +
										Old.get<V>(my) +
										Old.get<V>(py) +
										Old.get<V>(mx) +
                                        Old.get<V>(px) -
										6*Old.get<V>(Cp)) +
										deltaT * Old.get<U>(Cp) * Old.get<V>(Cp) * Old.get<V>(Cp) +
										- deltaT * (F+K) * Old.get<V>(Cp);

			// Next point in the grid
			++it;
		}

		//! \cond [stencil get and use] \endcond

		// Here we copy New into the old grid in preparation of the new step
		// It would be better to alternate, but using this we can show the usage
		// of the function copy. To note that copy work only on two grid of the same
		// decomposition. If you want to copy also the decomposition, or force to be
		// exactly the same, use Old = New
		Old.copy_sparse(New);

		// After copy we synchronize again the ghost part U and V
		Old.ghost_get<U,V>();

		// Every 500 time step we output the configuration for
		// visualization
		if (i % 500 == 0)
		{
			Old.save("output_" + std::to_string(count));
			count++;
		}
	}
	
	tot_sim.stop();
	std::cout << "Total simulation: " << tot_sim.getwct() << std::endl;

	//! \cond [time stepping] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse Gray Scott in 3D using sparse grids on CPU
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
	 * \page Grid_3_gs_3D_sparse Gray Scott in 3D using sparse grids on CPU
	 *
	 * # Full code # {#code}
	 *
	 * \include SparseGrid/1_gray_scott_3d_sparse/main.cpp
	 *
	 */
}
