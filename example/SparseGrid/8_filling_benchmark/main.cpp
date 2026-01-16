
#include "util/cuda_launch.hpp"
#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

/*!
 *
 * \page Grid_3_gs_3D_sparse_opt Gray Scott in 3D using sparse grids optimized on CPU
 *
 * [TOC]
 *
 * # Solving a gray scott-system in 3D using sparse grids optimized on CPU # {#1_gray_scott_3d_sparse}
 *
 * This example show how to solve a Gray-Scott system in 3D using sparse grids in an optimized way
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
 * This example is the same as \ref e3_gs_gray_scott_sparse the difference is optimizing for speed.
 *
 * Two optimization has been done. The first is to change the layout to struct of arrays defining the grid with
 *
 * \snippet SparseGrid/2_gray_scott_3d_sparse_opt/main.cpp grid definition
 *
 * The second is using the function **conv_cross2** to calculate the right-hand-side
 * this function can be used to do a convolution that involve points in a cross stencil like in figure that involve
 * two properties
 *
\verbatim

     *
     *
 * * x * *
     *
     *

\endverbatim
 *
 * The function accept a lambda function where the first 2 arguments are the output in form of Vc::double_v. If we use float
 * we have to use Vc::float_v or Vc::int_v in case the property is an integer. Vc variables come from the Vc library that is
 * now integrated in openfpm.
 *
 *\htmlonly
 * <a href="https://github.com/VcDevel/Vc" >Vc Library</a>
 *\endhtmlonly
 *
 * Vc::double_v in general pack 1,2,4 doubles dependently from the fact we choose to activate no-SSE,SSE or AVX at compiler level.
 * The arguments 3 and 4 contain the properties of two selected properties in the cross pattern given by xm xp ym yp zm zp.
 * The last arguments is instead the mask. The mask can be accessed to check the number of existing points. For example if
 * we have a cross stencil in 3D with stencil size = 1 than we expect 6 points. Note that the mask is an array because if Vc::double_v
 * contain 4 doubles than the mask has 4 elements accessed with the array operator []. The call **cross_conv2** also accept
 * template parameters the first two indicate the source porperties, the other two are the destination properties. While the
 * last is the extension of the stencil. In this case we use 1.
 *
 * The lambda function is defined as
 *
 * \snippet SparseGrid/2_gray_scott_3d_sparse_opt/main.cpp lambda
 *
 * and used in the body loop
 *
 * \snippet SparseGrid/2_gray_scott_3d_sparse_opt/main.cpp body
 *
 * To note that instead of copy we split the properties where we are acting at every iteration
 *
 */

constexpr int U = 0;
constexpr int V = 1;

constexpr int U_next = 2;
constexpr int V_next = 3;

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

void init(sgrid_dist_soa<3,double,aggregate<double,double,double,double> > & grid, Box<3,double> & domain)
{
	for (int i = 0 ; i < 10 ; i++)
	{
		timer t;
		t.start();

		auto it = grid.getGridIterator();

		while (it.isNext())
		{
			// Get the local grid key
			auto key = it.get_dist();

			// Old values U and V
			grid.template insert<U>(key) = 1.0;

			++it;
		}

		t.stop();
		std::cout << "Time populate: " << t.getwct() << std::endl;

		grid.clear();

                timer t2;
                t2.start();

                auto it2 = grid.getGridIterator();

                while (it2.isNext())
                {
                        // Get the local grid key
                        auto key = it2.get_dist();

                        // Old values U and V
                        grid.template insert<U>(key) = 5.0;
         
                        ++it2;
                }

		t2.stop();
                std::cout << "Time populate: " << t2.getwct() << std::endl;

	}
}


int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	// domain
	Box<3,double> domain({0.0,0.0,0.0},{2.5,2.5,2.5});
	
	// grid size
        size_t sz[3] = {256,256,256};

	// Define periodicity of the grid
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	
	// Ghost in grid unit
	Ghost<3,long int> g(1);
	
	// deltaT
	double deltaT = 0.25;

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

    //! \cond [grid definition] \endcond

	sgrid_dist_soa<3, double, aggregate<double,double,double,double>> grid(sz,domain,g,bc);

    //! \cond [grid definition] \endcond

	// spacing of the grid on x and y
	double spacing[3] = {grid.spacing(0),grid.spacing(1),grid.spacing(2)};

	init(grid,domain);

	// sync the ghost
	size_t count = 0;
	grid.template ghost_get<U,V>();

	// because we assume that spacing[x] == spacing[y] we use formula 2
	// and we calculate the prefactor of Eq 2
	double uFactor = deltaT * du/(spacing[x]*spacing[x]);
	double vFactor = deltaT * dv/(spacing[x]*spacing[x]);

	timer tot_sim;
	tot_sim.start();

	auto & v_cl = create_vcluster();
	
	tot_sim.stop();
	std::cout << "Total simulation: " << tot_sim.getwct() << std::endl;

	grid.write("final");

	//! \cond [time stepping] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse_opt Gray Scott in 3D using sparse grids optimized on CPU
	 *
	 * ## Finalize ##
	 *
	 * Deinitialize the library
	 *
	 * \snippet SparseGrid/2_gray_scott_3d_sparse_opt/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse_opt Gray Scott in 3D using sparse grids optimized on CPU
	 *
	 * # Full code # {#code}
	 *
	 * \include SparseGrid/2_gray_scott_3d_sparse_opt/main.cpp
	 *
	 */
}
