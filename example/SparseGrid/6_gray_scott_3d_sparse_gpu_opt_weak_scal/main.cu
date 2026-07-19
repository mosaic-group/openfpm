#include "Decomposition/Distribution/BoxDistribution.hpp"
#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

#ifdef CUDIFY_USE_METAL
using real_number = float;
#else
using real_number = double;
#endif

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

typedef sgrid_dist_id_gpu<3,float,aggregate<float,float,float,float>,CudaMemory, Dec> SparseGridType;

void init(SparseGridType & grid, Box<3,float> & domain, size_t (& div)[3])
{
	//! \cond [create points] \endcond

	real_number spacing_x = grid.spacing(0);
	real_number spacing_y = grid.spacing(1);
	real_number spacing_z = grid.spacing(2);

	typedef typename GetAddBlockType<SparseGridType>::type InsertBlockT;

	// Get the processor domain in continuos

	for (int i = 0 ; i < div[0] ; i++)
	{
		for (int j = 0 ; j < div[1] ; j++)
		{
			for (int k = 0 ; k < div[2] ; k++)
			{
				Point<3,real_number> p({0.5+i*1.0,0.5+j*1.0,0.5+k*1.0});
				Sphere<3,real_number> sph(p,0.3);

				Box<3,size_t> bx;

				for (int s = 0 ; s < 3 ; s++)
				{
					bx.setLow(s,(size_t)((sph.center(s) - 0.31)/grid.spacing(s)));
					bx.setHigh(s,(size_t)((sph.center(s) + 0.31)/grid.spacing(s)));
				}

				grid.addPoints([spacing_x,spacing_y,spacing_z,sph] __device__ (int i, int j, int k)
			        	{
							Point<3,real_number> pc({i*spacing_x,j*spacing_y,k*spacing_z});

							// Check if the point is in the domain
											if (sph.isInside(pc) )
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
				grid.removeUnusedBuffers();
			}
		}
	}

	for (int i = 0 ; i < div[0] ; i++)
	{
			for (int j = 0 ; j < div[1] ; j++)
			{
				Point<3,real_number> u({0.0,0.0,1.0});
				Point<3,real_number> c({0.5+i,0.5+j,0.0});

				Box<3,size_t> bx;

				bx.setLow(0,(0.4+i)/spacing_x);
				bx.setHigh(0,(0.6+i)/spacing_x);

				bx.setLow(1,(0.4+j)/spacing_y);
				bx.setHigh(1,(0.6+j)/spacing_y);

				bx.setLow(2,0);
				bx.setHigh(2,(size_t)grid.size(2));

				grid.addPoints(bx.getKP1(),bx.getKP2(),[spacing_x,spacing_y,spacing_z,u,c] __device__ (int i, int j, int k)
									{
						Point<3,real_number> pc({i*spacing_x,j*spacing_y,k*spacing_z});
												Point<3,real_number> pcs({i*spacing_x,j*spacing_y,k*spacing_z});
												Point<3,real_number> vp;

						// shift
						pc -= c; 

										// calculate the distance from the diagonal
										vp.get(0) = pc.get(1)*u.get(2) - pc.get(2)*u.get(1);
										vp.get(1) = pc.get(2)*u.get(0) - pc.get(0)*u.get(2);
										vp.get(2) = pc.get(0)*u.get(1) - pc.get(1)*u.get(0);

						real_number distance = vp.norm();

												// Check if the point is in the domain
												if (distance < 0.1 )
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
				grid.removeUnusedBuffers();
			}
	}

	for (int i = 0 ; i < div[0] ; i++)
	{
			for (int k = 0 ; k < div[2] ; k++)
			{
				Point<3,real_number> u({0.0,1.0,0.0});
				Point<3,real_number> c({0.5+i,0.0,0.5+k});

				Box<3,size_t> bx;

				bx.setLow(0,(0.4+i)/spacing_x);
				bx.setHigh(0,(0.6+i)/spacing_x);

				bx.setLow(2,(0.4+k)/spacing_z);
				bx.setHigh(2,(0.6+k)/spacing_z);

				bx.setLow(1,0);
				bx.setHigh(1,(size_t)grid.size(1));

				grid.addPoints(bx.getKP1(),bx.getKP2(),[spacing_x,spacing_y,spacing_z,u,c] __device__ (int i, int j, int k)
									{
						Point<3,real_number> pc({i*spacing_x,j*spacing_y,k*spacing_z});
												Point<3,real_number> pcs({i*spacing_x,j*spacing_y,k*spacing_z});
												Point<3,real_number> vp;

						// shift
						pc -= c; 

										// calculate the distance from the diagonal
										vp.get(0) = pc.get(1)*u.get(2) - pc.get(2)*u.get(1);
										vp.get(1) = pc.get(2)*u.get(0) - pc.get(0)*u.get(2);
										vp.get(2) = pc.get(0)*u.get(1) - pc.get(1)*u.get(0);

						real_number distance = vp.norm();

												// Check if the point is in the domain
												if (distance < 0.1 )
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
				grid.removeUnusedBuffers();
			}
	}

	for (int j = 0 ; j < div[1] ; j++)
	{
			for (int k = 0 ; k < div[2] ; k++)
			{
				Point<3,real_number> u({1.0,0.0,0.0});
				Point<3,real_number> c({0.0,0.5+j,0.5+k});

				Box<3,size_t> bx;

				bx.setLow(1,(0.4+j)/spacing_y);
				bx.setHigh(1,(0.6+j)/spacing_y);

				bx.setLow(2,(0.4+k)/spacing_z);
				bx.setHigh(2,(0.6+k)/spacing_z);

				bx.setLow(0,0);
				bx.setHigh(0,(size_t)grid.size(0));

				grid.addPoints(bx.getKP1(),bx.getKP2(),[spacing_x,spacing_y,spacing_z,u,c] __device__ (int i, int j, int k)
									{
						Point<3,real_number> pc({i*spacing_x,j*spacing_y,k*spacing_z});
												Point<3,real_number> pcs({i*spacing_x,j*spacing_y,k*spacing_z});
												Point<3,real_number> vp;

						// shift
						pc -= c; 

										// calculate the distance from the diagonal
										vp.get(0) = pc.get(1)*u.get(2) - pc.get(2)*u.get(1);
										vp.get(1) = pc.get(2)*u.get(0) - pc.get(0)*u.get(2);
										vp.get(2) = pc.get(0)*u.get(1) - pc.get(1)*u.get(0);

						real_number distance = vp.norm();

												// Check if the point is in the domain
												if (distance < 0.1 )
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
				grid.removeUnusedBuffers();
			}
	}

	//! \cond [create points] \endcond

	long int x_start = grid.size(0)*0.4f/domain.getHigh(0);
	long int y_start = grid.size(1)*0.4f/domain.getHigh(1);
	long int z_start = grid.size(1)*0.4f/domain.getHigh(2);

	long int x_stop = grid.size(0)*0.6f/domain.getHigh(0);
	long int y_stop = grid.size(1)*0.6f/domain.getHigh(1);
	long int z_stop = grid.size(1)*0.6f/domain.getHigh(2);

	//! \cond [create points sub] \endcond

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

	//! \cond [create points sub] \endcond
}


int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	// First we check which type of decomposition BoxDistritubion prodice
	auto & v_cl = create_vcluster();
	
	openfpm::vector<int> facts;
	getPrimeFactors(v_cl.size(),facts);

	size_t div[3];

	for (int i = 0 ; i < 3 ; i++)
	{div[i] = 1;}

	for (int i = 0 ; i < facts.size() ; i++)
	{div[i % 3] *= facts.get(i);}

	grid_sm<3,void> gdist(div);

	// domain
	Box<3,float> domain({0.0,0.0,0.0},{div[0]*1.0f,div[1]*1.0f,div[2]*1.0f});
	
	// grid size
	size_t sz[3] = {64*div[0],64*div[1],64*div[2]};

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

	SparseGridType grid(sz,domain,g,bc,0,gdist);

	// spacing of the grid on x and y
	float spacing[3] = {grid.spacing(0),grid.spacing(1),grid.spacing(2)};

	init(grid,domain,div);

        grid.deviceToHost<U,V,U_next,V_next>();
        grid.write("final");

	openfpm_finalize();
	return 0;

	// sync the ghost
	grid.template ghost_get<U,V>(RUN_ON_DEVICE);

	// because we assume that spacing[x] == spacing[y] we use formula 2
	// and we calculate the prefactor of Eq 2
	float uFactor = deltaT * du/(spacing[x]*spacing[x]);
	float vFactor = deltaT * dv/(spacing[x]*spacing[x]);

	timer tot_sim;
	tot_sim.start();

	for (size_t i = 0; i < timeSteps ; ++i)
	{
		if (v_cl.rank() == 0)
		{std::cout << "STEP: " << i << std::endl;}
/*		if (i % 300 == 0)
		{
			std::cout << "STEP: " << i << std::endl;
			grid.write_frame("out",i,VTK_WRITER);
		}*/

		//! \cond [stencil get and use] \endcond

		typedef typename GetCpBlockType<decltype(grid),0,1>::type CpBlockType;

		//! \cond [lambda] \endcond

		auto func = [uFactor,vFactor,deltaT,F,K] __device__ (float & u_out, float & v_out,
				                                   CpBlockType & u, CpBlockType & v,
				                                   int i, int j, int k){

				float uc = u(i,j,k);
				float vc = v(i,j,k);

				u_out = uc + uFactor *(u(i-1,j,k) + u(i+1,j,k) +
                                                       u(i,j-1,k) + u(i,j+1,k) +
                                                       u(i,j,k-1) + u(i,j,k+1) - 6.0f*uc) - deltaT * uc*vc*vc
                                                       - deltaT * F * (uc - 1.0f);


				v_out = vc + vFactor *(v(i-1,j,k) + v(i+1,j,k) +
                                                       v(i,j+1,k) + v(i,j-1,k) +
                                                       v(i,j,k-1) + v(i,j,k+1) - 6.0f*vc) + deltaT * uc*vc*vc
					               - deltaT * (F+K) * vc;
				};

		//! \cond [lambda] \endcond

		//! \cond [body] \endcond

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

		//! \cond [body] \endcond

		// Every 500 time step we output the configuration for
		// visualization
//		if (i % 500 == 0)
//		{
//			grid.save("output_" + std::to_string(count));
//			count++;
//		}
	}
	
	tot_sim.stop();
	std::cout << "Total simulation: " << tot_sim.getwct() << std::endl;

	grid.deviceToHost<U,V,U_next,V_next>();
	grid.write("final");

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
