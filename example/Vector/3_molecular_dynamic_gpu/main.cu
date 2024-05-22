/*!
 * \page Vector_3_md_dyn_gpu Vector 3 molecular dynamic on GPU
 *
 * [TOC]
 *
 * # Molecular Dynamic with Lennard-Jones potential on GPU {#e3_md_gpu}
 *
 * This Molecular dynamic simulation is exactly the same as \ref e3_md , with the difference that it work on Nvidia CUDA GPU and
 * for this reason we use one milions particles.
 * Before starting with this example, we suggest to start from the example \ref GPU_first_steps and \ref e3_md .
 *
 * \htmlonly
 * <img src="http://openfpm.mpi-cbg.de/web/images/examples/3_md_gpu/md_4_GPU.png"/>
 * \endhtmlonly
 *
 * # CUDA Kernels {#e3_cuda_ker_gpu}
 *
 * In the CPU example we have 4 loops
 *
 * * One loop to calculate forces
 * * One loop to do position and velocity integration as first step of the Verlet integration
 * * One loop to do velocity integration as second step of the Verlet integration
 * * One loop to calculate the total energy of the system
 *
 * All these loops must be converted into CUDA Kernels and launch the cuda kernel.
 *
 * ## Calculate forces {#e3_calc_force_gpu}
 *
 * The calculate forces function now has a call to the function \b getDomainIteratorGPU() \b, this function is equivalent
 *  to getDomainIterator(). What it return is struct with a convenient division in workgroups (\b wthr \b) and threads
 *  (\b thr \b) in a workgroups, to iterate across all the domain particles. The cuda kernel launch is a standard NVIDIA
 *  CUDA launch kernel, where if we want to pass the distributed data-structure to the kernel we have to remember to use the
 *  function toKernel().
 *
 * \snippet Vector/3_molecular_dynamic_gpu/main.cu calc_forces
 *
 * The kernel is templated on all the parameters passed with toKernel(). This is not strictly necessary but is convenient because
 * we let the compiler to infer the type returned by toKernel() without write the argument type explicitly. To get the particle index
 * we use the macro GET_PARTICLE this macro simply retrieve the particle id from the CUDA block id and thread is + it avoid that the
 * particle index does not overflow (The macro is simple to follow and is located in the file src/Vector/vector_dist_kernel.hpp).
 * The rest of the kernel is more or less a copy paste of the CPU iteration code
 *
 * \snippet Vector/3_molecular_dynamic_gpu/main.cu calculate_force_kernel
 *
 * For the two kernels that compose the 2 steps verlet integration, we simply followed the steps we described above. More interesting
 * is instead the calculation of the total energy of the system. In this case we have to do a full reduction, or accumulating the
 * energy of all the particles. In the CPU based code we used an external accumulator and we loop on all the particles summing
 * Kinetics energy and potential energy. On GPU this is not possible because such operation is serial and would not be fast on GPU.
 * To bypass this problem we have to first calculate the Kinetic energy + the potential energy for each particles in a kernel. This
 * kernel is very similar to the CPU code with the difference that we do not accumulate on an external variable, but we store
 * in a property per particle.
 *
 * \snippet Vector/3_molecular_dynamic_gpu/main.cu calculate_energy_kernel
 *
 * After that, we can use the function \b reduce_local<energy,\_add\_>(vd) \b to operate a parallel reduction on GPU of the property energy
 * on the vector_dist_gpu vd (the standard accumulation operation is the sum).
 *
 * \note reduce_local does only a local (within one GPU) reduction
 *
 * \snippet Vector/3_molecular_dynamic_gpu/main.cu calc_energy_red
 *
 * The rest remain mostly the same with the CPU code, with the exception that in the main loop where particle redistribution and ghost
 * fill happen directly on GPU we use the option \b RUN_ON_DEVICE \b.
 *
 * \snippet Vector/3_molecular_dynamic_gpu/main.cu run_on_device
 *
 * Every time we have to write a file we have to offload from GPU memory to host memory
 *
 * \snippet Vector/3_molecular_dynamic_gpu/main.cu move_to_host
 *
 */

#ifdef __NVCC__

#include "Vector/vector_dist.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include "timer.hpp"

#ifdef TEST_RUN
size_t nstep = 100;
#else
size_t nstep = 1000;
#endif

typedef float real_number;

//! \cond [constants] \endcond

constexpr int velocity = 0;
constexpr int force = 1;
constexpr int energy = 2;

//! \cond [calculate_force_kernel] \endcond

template<typename vector_dist_type,typename NN_type>
__global__ void calc_force_gpu(vector_dist_type vd, NN_type NN, real_number sigma12, real_number sigma6, real_number r_cut2)
{
	auto p = GET_PARTICLE(vd);

	// Get the position xp of the particle
	Point<3,real_number> xp = vd.getPos(p);

	// Reset the force counter
	vd.template getProp<force>(p)[0] = 0.0;
	vd.template getProp<force>(p)[1] = 0.0;
	vd.template getProp<force>(p)[2] = 0.0;



	// Get an iterator over the neighborhood particles of p
	auto Np = NN.getNNIterator(NN.getCell(vd.getPos(p)));

	// For each neighborhood particle ...
	while (Np.isNext())
	{
		// ... q
		auto q = Np.get();

		// if (p == q) skip this particle
		if (q == p)	{++Np; continue;};

		// Get the position of p
		Point<3,real_number> xq = vd.getPos(q);

		// Get the distance between p and q
		Point<3,real_number> r = xp - xq;

		// take the norm of this vector
		real_number rn = norm2(r);

		if (rn > r_cut2)
		{++Np; continue;};

		// Calculate the force, using pow is slower
		Point<3,real_number> f = 24.0*(2.0 *sigma12 / (rn*rn*rn*rn*rn*rn*rn) -  sigma6 / (rn*rn*rn*rn)) * r;

		// we sum the force produced by q on p
		vd.template getProp<force>(p)[0] += f.get(0);
		vd.template getProp<force>(p)[1] += f.get(1);
		vd.template getProp<force>(p)[2] += f.get(2);

		// Next neighborhood
		++Np;
	}
}

//! \cond [calculate_force_kernel] \endcond

template<typename vector_dist_type>
__global__ void update_velocity_position(vector_dist_type vd, real_number dt)
{
	auto p = GET_PARTICLE(vd);

	// here we calculate v(tn + 0.5)
	vd.template getProp<velocity>(p)[0] += 0.5*dt*vd.template getProp<force>(p)[0];
	vd.template getProp<velocity>(p)[1] += 0.5*dt*vd.template getProp<force>(p)[1];
	vd.template getProp<velocity>(p)[2] += 0.5*dt*vd.template getProp<force>(p)[2];

	// here we calculate x(tn + 1)
	vd.getPos(p)[0] += vd.template getProp<velocity>(p)[0]*dt;
	vd.getPos(p)[1] += vd.template getProp<velocity>(p)[1]*dt;
	vd.getPos(p)[2] += vd.template getProp<velocity>(p)[2]*dt;
}

template<typename vector_dist_type>
__global__ void update_velocity(vector_dist_type vd, real_number dt)
{
	auto p = GET_PARTICLE(vd);

	// here we calculate v(tn + 1)
	vd.template getProp<velocity>(p)[0] += 0.5*dt*vd.template getProp<force>(p)[0];
	vd.template getProp<velocity>(p)[1] += 0.5*dt*vd.template getProp<force>(p)[1];
	vd.template getProp<velocity>(p)[2] += 0.5*dt*vd.template getProp<force>(p)[2];
}

//! \cond [calculate_force_kernel] \endcond

//! \cond [calculate_energy_kernel] \endcond

template<typename vector_dist_type,typename NN_type>
__global__ void particle_energy(vector_dist_type vd, NN_type NN, real_number sigma12, real_number sigma6, real_number shift, real_number r_cut2)
{
	auto p = GET_PARTICLE(vd);

	// Get the position of the particle p
	Point<3,real_number> xp = vd.getPos(p);

	// Get an iterator over the neighborhood of the particle p
	auto Np = NN.getNNIterator(NN.getCell(vd.getPos(p)));

	real_number E = 0;

	// For each neighborhood of the particle p
	while (Np.isNext())
	{
		// Neighborhood particle q
		auto q = Np.get();

		// if p == q skip this particle
		if (q == p)	{++Np; continue;};

		// Get position of the particle q
		Point<3,real_number> xq = vd.getPos(q);

		// take the normalized direction
		real_number rn = norm2(xp - xq);

		if (rn > r_cut2)
		{++Np;continue;}

		// potential energy (using pow is slower)
		E += 2.0 * ( sigma12 / (rn*rn*rn*rn*rn*rn) - sigma6 / ( rn*rn*rn) ) - shift;

		// Next neighborhood
		++Np;
	}

	// Kinetic energy of the particle given by its actual speed
	vd.template getProp<energy>(p) = E + (vd.template getProp<velocity>(p)[0]*vd.template getProp<velocity>(p)[0] +
			vd.template getProp<velocity>(p)[1]*vd.template getProp<velocity>(p)[1] +
			vd.template getProp<velocity>(p)[2]*vd.template getProp<velocity>(p)[2]) / 2;
}

//! \cond [calculate_energy_kernel] \endcond

//! \cond [calc_forces] \endcond

template<typename CellList> void calc_forces(vector_dist_gpu<3,real_number, aggregate<real_number[3],real_number[3],real_number> > & vd, CellList & NN, real_number sigma12, real_number sigma6, real_number r_cut2)
{
	vd.updateCellListGPU(NN);

	// Get an iterator over particles
	auto it2 = vd.getDomainIteratorGPU();

	CUDA_LAUNCH(calc_force_gpu,it2,vd.toKernel(),NN.toKernel(),sigma12,sigma6,r_cut2);
}

//! \cond [calc_forces] \endcond

template<typename CellList> real_number calc_energy(vector_dist_gpu<3,real_number, aggregate<real_number[3],real_number[3],real_number> > & vd, CellList & NN, real_number sigma12, real_number sigma6, real_number r_cut2)
{
	real_number rc = r_cut2;
	real_number shift = 2.0 * ( sigma12 / (rc*rc*rc*rc*rc*rc) - sigma6 / ( rc*rc*rc) );

	vd.updateCellListGPU(NN);

	auto it2 = vd.getDomainIteratorGPU();

	CUDA_LAUNCH(particle_energy,it2,vd.toKernel(),NN.toKernel(),sigma12,sigma6,shift,r_cut2);

	//! \cond [calc_energy_red] \endcond

	// Calculated energy
	return reduce_local<energy,_add_>(vd);

	//! \cond [calc_energy_red] \endcond
}

int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	real_number sigma = 0.01;
	real_number r_cut = 3.0*sigma;

	// we will use it do place particles on a 10x10x10 Grid like
	size_t sz[3] = {100,100,100};

	// domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(r_cut);

	real_number dt = 0.00005;
	real_number sigma12 = pow(sigma,12);
	real_number sigma6 = pow(sigma,6);

	openfpm::vector<real_number> x;
	openfpm::vector<openfpm::vector<real_number>> y;

	vector_dist_gpu<3,real_number, aggregate<real_number[3],real_number[3],real_number> > vd(0,box,bc,ghost);

	// We create the grid iterator
	auto it = vd.getGridIterator(sz);

	while (it.isNext())
	{
		// Create a new particle
		vd.add();

		// key contain (i,j,k) index of the grid
		auto key = it.get();

		// The index of the grid can be accessed with key.get(0) == i, key.get(1) == j ...
		// We use getLastPos to set the position of the last particle added
		vd.getLastPos()[0] = key.get(0) * it.getSpacing(0);
		vd.getLastPos()[1] = key.get(1) * it.getSpacing(1);
		vd.getLastPos()[2] = key.get(2) * it.getSpacing(2);

		// We use getLastProp to set the property value of the last particle we added
		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;

		++it;
	}

	vd.hostToDevicePos();
	vd.hostToDeviceProp<velocity,force>();

	vd.map(RUN_ON_DEVICE);
	vd.ghost_get<>(RUN_ON_DEVICE);

	timer tsim;
	tsim.start();

	//! \cond [md steps] \endcond

	// Get the Cell list structure
	auto NN = vd.getCellListGPU(r_cut, CL_NON_SYMMETRIC, 2);

	// The standard
	// auto NN = vd.getCellList(r_cut);

	// calculate forces
	calc_forces(vd,NN,sigma12,sigma6,r_cut*r_cut);
	unsigned long int f = 0;

	// MD time stepping
	for (size_t i = 0; i < nstep ; i++)
	{
		// Get the iterator
		auto it3 = vd.getDomainIteratorGPU();

		CUDA_LAUNCH(update_velocity_position,it3,vd.toKernel(),dt);

		//! \cond [run_on_device] \endcond

		// Because we moved the particles in space we have to map them and re-sync the ghost
		vd.map(RUN_ON_DEVICE);
		vd.template ghost_get<>(RUN_ON_DEVICE);

		//! \cond [run_on_device] \endcond

		// calculate forces or a(tn + 1) Step 2
		calc_forces(vd,NN,sigma12,sigma6,r_cut*r_cut);


		// Integrate the velocity Step 3
		auto it4 = vd.getDomainIteratorGPU();

		CUDA_LAUNCH(update_velocity,it4,vd.toKernel(),dt);

		// After every iteration collect some statistic about the configuration
		if (i % 100 == 0)
		{
			//! \cond [move_to_host] \endcond

			vd.deviceToHostPos();
			vd.deviceToHostProp<0,1,2>();

			// We write the particle position for visualization (Without ghost)
			vd.deleteGhost();
			vd.write_frame("particles_",f);

			//! \cond [move_to_host] \endcond

			// we resync the ghost
			vd.ghost_get<>(RUN_ON_DEVICE);

			// We calculate the energy
			real_number energy = calc_energy(vd,NN,sigma12,sigma6,r_cut*r_cut);
			auto & vcl = create_vcluster();
			vcl.sum(energy);
			vcl.execute();

			// we save the energy calculated at time step i c contain the time-step y contain the energy
			x.add(i);
			y.add({energy});

			// We also print on terminal the value of the energy
			// only one processor (master) write on terminal
			if (vcl.getProcessUnitID() == 0)
				std::cout << "Energy: " << energy << std::endl;

			f++;
		}
	}

	//! \cond [md steps] \endcond

	tsim.stop();
	std::cout << "Time: " << tsim.getwct() << std::endl;

	// Google charts options, it store the options to draw the X Y graph
	GCoptions options;

	// Title of the graph
	options.title = std::string("Energy with time");

	// Y axis name
	options.yAxis = std::string("Energy");

	// X axis name
	options.xAxis = std::string("iteration");

	// width of the line
	options.lineWidth = 1.0;

	// Resolution in x
	options.width = 1280;

	// Resolution in y
	options.heigh = 720;

	// Add zoom capability
	options.more = GC_ZOOM;

	// Object that draw the X Y graph
	GoogleChart cg;

	// Add the graph
	// The graph that it produce is in svg format that can be opened on browser
	cg.AddLinesGraph(x,y,options);

	// Write into html format
	cg.write("gc_plot2_out.html");

	openfpm_finalize();
}

#else

int main(int argc, char* argv[])
{
        return 0;
}

#endif


