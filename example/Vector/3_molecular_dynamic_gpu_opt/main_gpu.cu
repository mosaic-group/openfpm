/*!
 * \page Vector_3_md_dyn_gpu_opt Vector 3 molecular dynamic on GPU (optimized version)
 *
 * [TOC]
 *
 * # Molecular Dynamic with Lennard-Jones potential on GPU (Optimized) {#e3_md_gpu_opt}
 *
 * \htmlonly
 * <img src="http://openfpm.mpi-cbg.de/web/images/examples/3_md_gpu/md_4_GPU.png"/>
 * \endhtmlonly
 *
 * In this optimized version of the example \ref e3_md_gpu, we operate two optimization:
 *
 * * We make the access coalesced
 * * We use half radius cell-list spacing
 *
 * ## Coalesced access {#e3_gpu_opt_ca}
 *
 * In GPU to get the maximum performance it is very important to access in a coalesced way. Access in a coalesced way mean that
 *  if thread 1 access adress 0x000 thread 2 (in the same Streaming Multiprocessors) should ideally access 0x004 or more in
 *  general an adress in the same cache line. Another factor that contribute to speed is to overall restrict the threads in
 *  the same SM should to possibly work on a limited number of caches lines so that the L1 cache of each SM could optimally
 *  speed up the access to global memory.
 *
 *  Unfortunately particles by nature can be randomly distributed in space and memory, and reach the ideal situation
 * in the case of neighborhood access of the particles can be challenging. Suppose that thread 1 take particle 0 and thread 2
 * take particle 1, but 0 and 1 are far in space the neighborhood of 0 does not overlap the neighborhood of 1. This mean that
 * most probably the access of the neighborhood of 0 will be scattered in memory and the same for 1, having an extremely low
 * probability that 2 thread in the SM hit the same cache line and increasing the number of cache lines the SM has to retrieve.
 *  On the other hand if we reorder the particle in memory by their spatial position using a cell list like in figure.
 *
 * \htmlonly
  <table style="width:100%">
  <tr>
    <td><img src="http://openfpm.mpi-cbg.de/web/images/examples/3_md_gpu/vector_sorted.jpg"/></td>
    <td><img src="http://openfpm.mpi-cbg.de/web/images/examples/3_md_gpu/vector_unsorted.jpg"/></td>
  </tr>
  <tr>
    <td><strong>Fig1: Sorted vector</strong></td>
    <td><strong>Fig2: Unsorted vector</strong></td>
  </tr>
</table>
  \endhtmlonly
 *
 * We can see that now the neighborhood of particle 0 and particle 1 overlap increasing that chance of cache hit, additionally if
 * all particles processed by one SM stay in one cell or few neighborhood cell, the number of cache line that an SM has to read is
 * reduced, with a significant speed-up.
 *
 * In OpenFPM \b getCellListGPU() \b produces a re-ordered version of the original vector when CL_GPU_REORDER option is enabled.
 *
 * \snippet Vector/3_molecular_dynamic_gpu_opt/main_gpu.cu calc_force_sorted
 *
 * The rest remain mainly the same, with the expectation, that we now use the macro GET_PARTICLE_SORT. This macro is similar to GET_PARTICLE
 * but with a substantial difference. While in the normal unsorted vector particles in the ghost area are always added at the end
 * in the sorted one domain + ghost are reordered, and there is not a clear separation between them. This mean that we need a list of all
 * the domain particles, if we want iterate cross them. GET_PARTICLE_SORT use a list to convert thread index to domain particle index.
 * Additionally when we get a neighborhood iterator from the Cell-list we must use \bget_sorted_index\b instead of \bget\b
 *
 * \snippet Vector/3_molecular_dynamic_gpu_opt/main_gpu.cu get_sorted_index
 *
 * After we launched the kernel all the data are written in the sorted vector. In order to merge back the data to the unsorted one
 * we have to use the function \b vd.merge_sort<force>(NN) \b. Where vd is the vector_dist_gpu where we want to merge the
 * data from sorted to non sorted. \b force \b is the property we want to merge and \b NN \b is the Cell-list that produced the
 * sorted distribution.
 *
 * \snippet Vector/3_molecular_dynamic_gpu_opt/main_gpu.cu merge_sort
 *
 * \note it is possible to launch multiple kernel on the sorted version, but consider that at some point the data must be merged
 * back because functions like map and ghost_get work on the unsorted version
 *
 * ## Half radius cell-list spacing {#e3_gpu_opt_hr}
 *
 * Using Cell-lists with spacing equal to the radius in general require to fetch all the 9 cells in 2D and 27 cells in 3D. All the
 * particles in such cells include particles within radius r and others more distant than r. This mean that we have to filter the particles
 * checking the radius. It is possible to filter further more the particles using finer cell-list cells. Suppose that you use
 * cell-lists with spacing half of the radius. we just the to check the 25 cells in 2D and the 125 cells in 3D. While we have more
 * cells the overall volume spanned by the 25/125 cells is just a fraction. In fact the surface of the 25 cells is given by
 *
 * \f$ (5\frac{h}{2})^2 = \frac{25}{4} h^2 \f$
 * \f$ (5\frac{h}{2})^3 = \frac{125}{8} h^3 \f$
 *
 * while for the normal cell-list is
 *
 * \f$ (3h)^2 = 9h^2 \f$
 * \f$ (3h)^3 = 27h^3 \f$
 *
 * This mean that the finer cell-list in order to find the neighborhood particles use an area smaller: precisely is 69% of
 * the normal cell-list in 2D, and 57% of the normal cell-list in 3D. In particles this mean that normal cell-list return
 * in average 45% more particles in 2D and  75% more in 3D.
 *
 * Constructing an half spacing cell-list is standard. In the function \b getCellListGPU \b we specify half radius
 *
 * \snippet Vector/3_molecular_dynamic_gpu_opt/main_gpu.cu get_half_cl
 *
 * while to use it, instead of the \b getNNIteratorBox \b we use
 *
 * \note \b getNNIteratorBox \b is initialized in cell-list by \b setBoxNN \b (default = 1) before it transferred to gpu that indicate how many neighborhood cell the NN iterator has to span.
 *       For example \b setBoxNN(1) \b is the standard 9/27 neighborhood cell-list.\b setBoxNN(2) \b is the 25/125 neighborhood
 *        and so on.
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

constexpr int velocity = 0;
constexpr int force = 1;
constexpr int energy = 2;

template<typename vector_dist_type,typename NN_type>
__global__ void calc_force_gpu(vector_dist_type vd, NN_type NN, real_number sigma12, real_number sigma6, real_number r_cut2)
{
	unsigned int p;
	GET_PARTICLE_SORT(p,NN);

	// Get the position xp of the particle
	Point<3,real_number> xp = vd.getPos(p);

	// Reset the force counter
	vd.template getProp<force>(p)[0] = 0.0;
	vd.template getProp<force>(p)[1] = 0.0;
	vd.template getProp<force>(p)[2] = 0.0;

	Point<3,real_number> force_;
	force_.get(0) = 0.0;
	force_.get(1) = 0.0;
	force_.get(2) = 0.0;

	// Get an iterator over the neighborhood particles of p
	auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(p)));

	// For each neighborhood particle ...
	while (Np.isNext())
	{
		//! \cond [get_sorted_index] \endcond

		// ... q
		auto q = Np.get_sort();

		//! \cond [get_sorted_index] \endcond

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
		force_ += f;

		// Next neighborhood
		++Np;
	}

	// we sum the force produced by q on p
	vd.template getProp<force>(p)[0] = force_.get(0);
	vd.template getProp<force>(p)[1] = force_.get(1);
	vd.template getProp<force>(p)[2] = force_.get(2);
}

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

template<typename vector_dist_type,typename NN_type>
__global__ void particle_energy(vector_dist_type vd, NN_type NN, real_number sigma12, real_number sigma6, real_number shift, real_number r_cut2)
{
	unsigned int p;
	GET_PARTICLE_SORT(p,NN);

	// Get the position of the particle p
	Point<3,real_number> xp = vd.getPos(p);

	// Get an iterator over the neighborhood of the particle p
	auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(p)));

	real_number E = 0;

	// For each neighborhood of the particle p
	while (Np.isNext())
	{
		// Neighborhood particle q
		auto q = Np.get_sort();

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

template<typename CellList> void calc_forces(vector_dist_gpu<3,real_number, aggregate<real_number[3],real_number[3],real_number> > & vd, CellList & NN, real_number sigma12, real_number sigma6, real_number r_cut2)
{
	// Get an iterator over particles
	auto it2 = vd.getDomainIteratorGPU();

	//! \cond [calc_force_sorted] \endcond

	// reorder positions only (no properties)
	// as nothing else is needed to be read in calc_force_gpu
	vd.template updateCellListGPU<>(NN);

	CUDA_LAUNCH(calc_force_gpu,it2,vd.toKernel(),NN.toKernel(),sigma12,sigma6,r_cut2);

	//! \cond [calc_force_sorted] \endcond
	vd.template restoreOrder<force>(NN);
}

template<typename CellList> real_number calc_energy(vector_dist_gpu<3,real_number, aggregate<real_number[3],real_number[3],real_number> > & vd, CellList & NN, real_number sigma12, real_number sigma6, real_number r_cut2)
{
	real_number rc = r_cut2;
	real_number shift = 2.0 * ( sigma12 / (rc*rc*rc*rc*rc*rc) - sigma6 / ( rc*rc*rc) );

	vd.template updateCellListGPU<velocity>(NN);

	auto it2 = vd.getDomainIteratorGPU();

	CUDA_LAUNCH(particle_energy,it2,vd.toKernel(),NN.toKernel(),sigma12,sigma6,shift,r_cut2);

	vd.template restoreOrder<energy>(NN);

	// Calculated energy
	return reduce_local<energy,_add_>(vd);
}

int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	real_number sigma = 0.01;
	real_number r_cut =3.0*sigma;

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

	//! \cond [get_half_cl] \endcond

	// Get the Cell list structure
	auto NN = vd.getCellListGPU(r_cut / 2.0, CL_NON_SYMMETRIC | CL_GPU_REORDER, 2);

	//! \cond [get_half_cl] \endcond

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

		// Because we moved the particles in space we have to map them and re-sync the ghost
		vd.map(RUN_ON_DEVICE);
		vd.template ghost_get<>(RUN_ON_DEVICE);

		// calculate forces or a(tn + 1) Step 2
		calc_forces(vd,NN,sigma12,sigma6,r_cut*r_cut);

		// Integrate the velocity Step 3
		auto it4 = vd.getDomainIteratorGPU();

		CUDA_LAUNCH(update_velocity,it4,vd.toKernel(),dt);

		// After every iteration collect some statistic about the configuration
		if (i % 100 == 0)
		{
			vd.deviceToHostPos();
			vd.deviceToHostProp<0,1,2>();

			// We write the particle position for visualization (Without ghost)
			vd.deleteGhost();
			vd.write_frame("particles_",f);

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


