#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

/*!
 *
 * \page Grid_3_gs_3D_sparse_gpu_cs Gray Scott in 3D using sparse grids on GPU in complex geometry
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

typedef sgrid_dist_id_gpu<3,double,aggregate<double,double,double,double> > sgrid_type;

void init(sgrid_type & grid, Box<3,double> & domain)
{
	auto it = grid.getGridIterator();
	Point<3,double> p1;
	Point<3,double> p2;
	Point<3,double> p3;

	// Shere1
	for (int i = 0 ; i < 3 ; i++)
	{p1.get(i) = 2.0;}

        // Shere2
        for (int i = 0 ; i < 3 ; i++)
        {p2.get(i) = 1.0;}

        // Shere3
        for (int i = 0 ; i < 3 ; i++)
        {p3.get(i) = 0.5;}

	Sphere<3,double> sph1(p1,0.3);
	Sphere<3,double> sph2(p2,0.3);
	Sphere<3,double> sph3(p3,0.3);

	Point<3,double> u({1.0,1.0,1.0});
	Box<3,double> channel_box(p3,p1);

	double spacing_x = grid.spacing(0);
	double spacing_y = grid.spacing(1);
	double spacing_z = grid.spacing(2);

	typedef typename GetAddBlockType<sgrid_type>::type InsertBlockT;

	grid.addPoints([spacing_x,spacing_y,spacing_z,u,sph1,sph2,sph3,channel_box] __device__ (int i, int j, int k)
			        {
						Point<3,double> pc({i*spacing_x,j*spacing_y,k*spacing_z});
						Point<3,double> vp;

		                // calculate the distance from the diagonal
		                vp.get(0) = pc.get(1)*u.get(2) - pc.get(2)*u.get(1);
		                vp.get(1) = pc.get(2)*u.get(0) - pc.get(0)*u.get(2);
		                vp.get(2) = pc.get(0)*u.get(1) - pc.get(1)*u.get(0);

		                double distance = vp.norm() / sqrt(3.0f);

		                // Check if the point is in the domain
		                if (sph1.isInside(pc) || sph2.isInside(pc) || sph3.isInside(pc) || (distance < 0.1 && channel_box.isInside(pc)) )
		                {return true;}

						return false;
			        },
			        [] __device__ (InsertBlockT & data, int i, int j, int k)
			        {
			        	data.template get<U>() = 1.0;
			        	data.template get<V>() = 0.0;
			        }
			        );


	grid.template flush<smax_<U>,smax_<V>>(flush_type::FLUSH_ON_DEVICE);

	long int x_start = grid.size(0)*1.95f/domain.getHigh(0);
	long int y_start = grid.size(1)*1.95f/domain.getHigh(1);
	long int z_start = grid.size(1)*1.95f/domain.getHigh(2);

	long int x_stop = grid.size(0)*2.05f/domain.getHigh(0);
	long int y_stop = grid.size(1)*2.05f/domain.getHigh(1);
	long int z_stop = grid.size(1)*2.05f/domain.getHigh(2);

	grid_key_dx<3> start({x_start,y_start,z_start});
	grid_key_dx<3> stop ({x_stop,y_stop,z_stop});

        grid.addPoints(start,stop,[] __device__ (int i, int j, int k)
                                {
                                                return true;
                                },
                                [] __device__ (InsertBlockT & data, int i, int j, int k)
                                {
                                        data.template get<U>() = 0.5;
                                        data.template get<V>() = 0.24;
                                }
                                );

	grid.template flush<smin_<U>,smax_<V>>(flush_type::FLUSH_ON_DEVICE);
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
	double deltaT = 0.1;

	// Diffusion constant for specie U
	double du = 2*1e-5;

	// Diffusion constant for specie V
	double dv = 1*1e-5;

#ifdef TEST_RUN
        // Number of timesteps
        size_t timeSteps = 300;
#else
	// Number of timesteps
        size_t timeSteps = 150000;
#endif

	// K and F (Physical constant in the equation)
        double K = 0.053;
        double F = 0.014;

	sgrid_type grid(sz,domain,g,bc);

	grid.template setBackgroundValue<0>(-0.5);
	grid.template setBackgroundValue<1>(-0.5);
	grid.template setBackgroundValue<2>(-0.5);
	grid.template setBackgroundValue<3>(-0.5);
	
	// spacing of the grid on x and y
	double spacing[3] = {grid.spacing(0),grid.spacing(1),grid.spacing(2)};

	init(grid,domain);

	// sync the ghost
	grid.template ghost_get<U,V>(RUN_ON_DEVICE);

	// because we assume that spacing[x] == spacing[y] we use formula 2
	// and we calculate the prefactor of Eq 2
	double uFactor = deltaT * du/(spacing[0]*spacing[0]);
	double vFactor = deltaT * dv/(spacing[0]*spacing[0]);

	grid.template deviceToHost<U,V>();
	grid.write("Init_condition");

	timer tot_sim;
	tot_sim.start();

	for (size_t i = 0; i < timeSteps; ++i)
	{
		//! \cond [stencil get and use] \endcond

        		typedef typename GetCpBlockType<decltype(grid),0,1>::type CpBlockType;

        		auto func = [uFactor,vFactor,deltaT,F,K] __device__ (double & u_out, double & v_out,
        				                                   CpBlockType & u, CpBlockType & v,
        				                                   int i, int j, int k){

        				double uc = u(i,j,k);
        				double vc = v(i,j,k);

        				double u_px = u(i+1,j,k);
        				double u_mx = u(i-1,j,k);

        				double u_py = u(i,j+1,k);
        				double u_my = u(i,j-1,k);

        				double u_pz = u(i,j,k+1);
        				double u_mz = u(i,j,k-1);

        				double v_px = v(i+1,j,k);
        				double v_mx = v(i-1,j,k);

        				double v_py = v(i,j+1,k);
        				double v_my = v(i,j-1,k);

        				double v_pz = v(i,j,k+1);
        				double v_mz = v(i,j,k-1);

        				// U fix

        				if (u_mx < -0.1 && u_px < -0.1)
        				{
        					u_mx = uc;
        					u_px = uc;
        				}

        				if (u_mx < -0.1)
        				{u_mx = u_px;}

        				if (u_px < -0.1)
        				{u_px = u_mx;}

        				if (u_my < -0.1 && u_py < -0.1)
        				{
        					u_my = uc;
        					u_py = uc;
        				}

        				if (u_my < -0.1)
        				{u_my = u_py;}

        				if (u_py < -0.1)
        				{u_py = u_my;}

        				if (u_mz < -0.1 && u_pz < -0.1)
        				{
        					u_mz = uc;
        					u_pz = uc;
        				}

        				if (u_mz < -0.1)
        				{u_mz = u_pz;}

        				if (u_pz < -0.1)
        				{u_pz = u_mz;}

        				// V fix

        				if (v_mx < -0.1 && v_px < -0.1)
        				{
        					v_mx = uc;
        					v_px = uc;
        				}

        				if (v_mx < -0.1)
        				{v_mx = v_px;}

        				if (v_px < -0.1)
        				{v_px = v_mx;}

        				if (v_my < -0.1 && v_py < -0.1)
        				{
        					v_my = uc;
        					v_py = uc;
        				}

        				if (v_my < -0.1)
        				{v_my = v_py;}

        				if (v_py < -0.1)
        				{v_py = v_my;}

        				if (v_mz < -0.1 && v_pz < -0.1)
        				{
        					v_mz = uc;
        					v_pz = uc;
        				}

        				if (v_mz < -0.1)
        				{v_mz = v_pz;}

        				if (v_pz < -0.1)
        				{v_pz = v_mz;}

        				u_out = uc + uFactor *(u_mx + u_px +
                                                               u_my + u_py +
                                                               u_mz + u_pz - 6.0*uc) - deltaT * uc*vc*vc
                                                               - deltaT * F * (uc - 1.0);


        				v_out = vc + vFactor *(v_mx + v_px +
                                                               v_py + v_my +
                                                               v_mz + v_pz - 6.0*vc) + deltaT * uc*vc*vc
        					               - deltaT * (F+K) * vc;

        				};


        		if (i % 2 == 0)
        		{
        			grid.conv2<U,V,U_next,V_next,1>({0,0,0},{(long int)sz[0]-1,(long int)sz[1]-1,(long int)sz[2]-1},func);

        			// After copy we synchronize again the ghost part U and V

        			grid.ghost_get<U_next,V_next>(RUN_ON_DEVICE | SKIP_LABELLING);
        		}
        		else
        		{
        			grid.conv2<U_next,V_next,U,V,1>({0,0,0},{(long int)sz[0]-1,(long int)sz[1]-1,(long int)sz[2]-1},func);

        			// After copy we synchronize again the ghost part U and V
        			grid.ghost_get<U,V>(RUN_ON_DEVICE | SKIP_LABELLING);
        		}

		//! \cond [stencil get and use] \endcond

		// After copy we synchronize again the ghost part U and V

		// Every 500 time step we output the configuration for
		// visualization
/*		if (i % 500 == 0)
		{
			grid.save("output_" + std::to_string(count));
			count++;
		}*/

                std::cout << "STEP: " << i  << std::endl;
/*                if (i % 300 == 0)
                {
                	grid.template deviceToHost<U,V>();
                        grid.write_frame("out",i);
                }*/
	}
	
	tot_sim.stop();
	std::cout << "Total simulation: " << tot_sim.getwct() << std::endl;

	grid.template deviceToHost<U,V>();
	grid.write("Final");

	//! \cond [time stepping] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse_gpu_cs Gray Scott in 3D using sparse grids on GPU in complex geometry
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
	 * \page Grid_3_gs_3D_sparse_gpu_cs Gray Scott in 3D using sparse grids on GPU in complex geometry
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

