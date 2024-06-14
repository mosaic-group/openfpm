
#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include "timer.hpp"

/*!
 * \page Vector_4_reo_root Vector 4 reordering
 *   \subpage Vector_4_reo
 *   \subpage Vector_4_comp_reo
 *
 */

/*!
 * \page Vector_4_comp_reo Vector 4 computational reordering and cache friendliness
 *
 * [TOC]
 *
 * # Computation-reordering and cache friendly computation # {#e4_comp_reo}
 *
 * In this example we show how reordering the data can significantly improve the computation speed.
 * In order to do this we will re-work the molecular dynamic example.
 *
 *
 * ## Calculate forces ##
 *
 * This function is the same as the molecular dynamic example with few changes:
 *
 * * The function now take as argument CellList_hilb instead of CellList
 *   \snippet Vector/4_reorder/main_comp_ord.cpp calc forces
 *
 * * We get an iterator from the Cell list instead that from the vector
 *   \snippet Vector/4_reorder/main_comp_ord.cpp NN iterator
 *
 * \see \ref e3_md_cf
 *
 */

//! \cond [calc forces all] \endcond

constexpr int velocity = 0;
constexpr int force = 1;

//! \cond [calc forces] \endcond

template<typename CellList> void calc_forces(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, CellList & NN, double sigma12, double sigma6)
{

//! \cond [calc forces] \endcond


	// Uodate the cell-list
	vd.updateCellList(NN);

	//! \cond [NN iterator] \endcond

	// Get an iterator over particles
	auto it2 = NN.getCellParticleIterator();

	//! \cond [NN iterator] \endcond

	// For each particle p ...
	while (it2.isNext())
	{
		// ... get the particle p
		auto p = it2.get();

		// Get the position xp of the particle
		Point<3,double> xp = vd.getPos(p);

		// Reset the force counter
		vd.template getProp<force>(p)[0] = 0.0;
		vd.template getProp<force>(p)[1] = 0.0;
		vd.template getProp<force>(p)[2] = 0.0;

		// Get an iterator over the neighborhood particles of p
		auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(p)));

		// For each neighborhood particle ...
		while (Np.isNext())
		{
			// ... q
			auto q = Np.get();

			// if (p == q) skip this particle
			if (q == p)	{++Np; continue;};

			// Get the position of p
			Point<3,double> xq = vd.getPos(q);

			// Get the distance between p and q
			Point<3,double> r = xp - xq;

			// take the norm of this vector
			double rn = norm2(r);

			// Calculate the force, using pow is slower
			Point<3,double> f = 24.0*(2.0 *sigma12 / (rn*rn*rn*rn*rn*rn*rn) -  sigma6 / (rn*rn*rn*rn)) * r;

			// we sum the force produced by q on p
			vd.template getProp<force>(p)[0] += f.get(0);
			vd.template getProp<force>(p)[1] += f.get(1);
			vd.template getProp<force>(p)[2] += f.get(2);

			// Next neighborhood
			++Np;
		}

		// Next particle
		++it2;
	}
}

//! \cond [calc forces all] \endcond

/*!
 *
 * \page Vector_4_comp_reo Vector 4 computational reordering and cache friendliness
 *
 *
 *
 * ## Calculate energy ##
 *
 * In this function is the same as the molecular dynamic example with few changes:
 *
 * * The function now take as argument CellList_hilb instead of CellList
 *   \snippet Vector/4_reorder/main_comp_ord.cpp calc forces
 *
 * * We get an iterator from the Cell list instead that from the vector
 *   \snippet Vector/4_reorder/main_comp_ord.cpp NN iterator
 *   The difference in doing this is that now we iterate on particles in a smarter way.
 *   We will explain more in detail later in the example
 *
 *
 * \see \ref e3_md_ce
 *
 */

//! \cond [calc energy all] \endcond

template<typename CellList> double calc_energy(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, CellList & NN, double sigma12, double sigma6)
{
	double E = 0.0;

	// update cell-list
	vd.updateCellList(NN);

	// Get the iterator
	auto it2 = NN.getCellParticleIterator();

	// For each particle ...
	while (it2.isNext())
	{
		// ... p
		auto p = it2.get();

		// Get the position of the particle p
		Point<3,double> xp = vd.getPos(p);

		// Reset the force
		vd.template getProp<force>(p)[0] = 0.0;
		vd.template getProp<force>(p)[1] = 0.0;
		vd.template getProp<force>(p)[2] = 0.0;

		// Get an iterator over the neighborhood of the particle p
		auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(p)));

		// For each neighborhood of the particle p
		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			// if p == q skip this particle
			if (q == p)	{++Np; continue;};

			// Get position of the particle q
			Point<3,double> xq = vd.getPos(q);

			// take the normalized direction
			double rn = norm2(xp - xq);

			// potential energy (using pow is slower)
			E += 4.0 * ( sigma12 / (rn*rn*rn*rn*rn*rn) - sigma6 / ( rn*rn*rn) );

			// Next neighborhood
			++Np;
		}

		// Kinetic energy of the particle given by its actual speed
		E +=   (vd.template getProp<velocity>(p)[0]*vd.template getProp<velocity>(p)[0] +
				vd.template getProp<velocity>(p)[1]*vd.template getProp<velocity>(p)[1] +
				vd.template getProp<velocity>(p)[2]*vd.template getProp<velocity>(p)[2]) / 2;

		// Next Particle
		++it2;
	}

	// Calculated energy
	return E;
}

//! \cond [calc energy all] \endcond


int main(int argc, char* argv[])
{
	/*!
	 *
	 * \page Vector_4_comp_reo Vector 4 computational reordering and cache friendliness
	 *
	 * ## Initialization ##
	 *
	 * The initialization is the same as the molecular dynamic example. The differences are in the
	 * parameters. We will use a bigger system, with more particles. The delta time for integration
	 * is chosen in order to keep the system stable.
	 *
	 * \see \ref e3_md_init
	 *
	 * \snippet Vector/4_reorder/main_comp_ord.cpp vect create
	 *
	 */

	//! \cond [vect create] \endcond

	double dt = 0.0001;
	double r_cut = 0.03;
	double sigma = r_cut/3.0;
	double sigma12 = pow(sigma,12);
	double sigma6 = pow(sigma,6);

	openfpm::vector<double> x;
	openfpm::vector<openfpm::vector<double>> y;

	openfpm_init(&argc,&argv);
	Vcluster<> & v_cl = create_vcluster();

	// we will use it do place particles on a 40x40x40 Grid like
	size_t sz[3] = {40,40,40};

	// domain
	Box<3,double> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,double> ghost(r_cut);

	vector_dist<3,double, aggregate<double[3],double[3]> > vd(0,box,bc,ghost);

	//! \cond [vect create] \endcond

	/*!
	 * \page Vector_4_comp_reo Vector 4 computational reordering and cache friendliness
	 *
	 * ## Particles on a grid like position ##
	 *
	 * Here we place the particles on a grid like manner
	 *
	 * \see \ref e3_md_gl
	 *
	 * \snippet Vector/4_reorder/main_comp_ord.cpp vect grid
	 *
	 */

	//! \cond [vect grid] \endcond

	auto it = vd.getGridIterator(sz);

	while (it.isNext())
	{
		vd.add();

		auto key = it.get();

		vd.getLastPos()[0] = key.get(0) * it.getSpacing(0);
		vd.getLastPos()[1] = key.get(1) * it.getSpacing(1);
		vd.getLastPos()[2] = key.get(2) * it.getSpacing(2);

		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;

		++it;
	}

	//! \cond [vect grid] \endcond

	/*!
	 *
	 * \page Vector_4_comp_reo Vector 4 computational reordering and cache friendliness
	 *
	 * ## Molecular dynamic steps ##
	 *
	 * Here we do 30000 MD steps using verlet integrator the cycle is the same as the
	 * molecular dynamic example. with the following changes.
	 *
	 * ### Cell lists ###
	 *
	 * Instead of getting the normal cell list we get an hilbert curve cell-list by passing flag **CL_HILBERT_CELL_KEYS**
	 * to **getCellList** or by directly calling **getCellList_hilb**. Such cell list has a
	 * function called **getCellParticleIterator** used inside the function **calc_forces** and **calc_energy**
	 * that iterate across all the particles but in a smart-way. In practice
	 * given an r-cut a cell-list is constructed with the provided spacing. Suppose to have a cell-list
	 * \f$ m \times n \f$, an hilbert curve \f$ 2^k \times 2^k \f$ is contructed with \f$ k = ceil(log_2(max(m,n))) \f$.
	 * Cell-lists are explored according to this Hilbert curve, If a cell does not exist is simply skipped.
	 *
	 *
	 * \verbatim
	+------+------+------+------+     Example of Hilbert curve running on a 3 x 3 Cell
	|      |      |      |      |     An hilbert curve of k = ceil(log_2(3)) = 4
	|  X+---->X   |  X +---> X  |
	|  ^   |  +   |  ^   |   +  |
	***|******|******|****---|--+      *******
	*  +   |  v   |  +   *   v  |      *     *
	*  7   |  8+---->9   *   X  |      *     *  = Domain
	*  ^   |      |      *   +  |      *     *
	*--|-----------------*---|--+      *******
	*  +   |      |      *   v  |
	*  4<----+5   |   6<---+ X  |
	*      |  ^   |   +  *      |
	*---------|-------|--*------+
	*      |  +   |   v  *      |
	*  1+---->2   |   3+---> X  |
	*      |      |      *      |
	**********************------+

     this mean that we will iterate the following cells

     1,2,5,4,7,8,9,6,3

	 Suppose now that the particles are ordered like described




    Particles   id      Cell
                 0         1
                 1         7
                 2         8
				 3		   1
				 4		   9
				 5		   9
				 6		   6
				 7		   7
				 8		   3
				 9		   2
				10		   4
				11		   3


	The iterator of the cell-list will explore the particles in the following way

	Cell     1  2 5 4  7  8  9  6 3
       	   |   | | | |   | |   | | |
	    	0,3,9,,10,1,7,2,4,5,6,8


	 * \endverbatim
	 *
	 * We cannot explain here what is a cache, but in practice is a fast memory in the CPU able
	 * to store chunks of memory. The cache in general is much smaller than RAM, but the big advantage
	 * is its speed. Retrieve data from the cache is much faster than RAM. Unfortunately the factors
	 *  that determine what is on cache and what is not are multiples: Type of cache, algorithm ... .
	 * Qualitatively all caches will tend to load chunks of data that you read multiple-time, or chunks
	 *  of data that probably you will read based on pattern analysis. A small example is a linear memory copy where
	 * you read consecutively memory and you write on consecutive memory.
	 * Modern CPU recognize such pattern and decide to load on cache the consecutive memory before
	 *  you actually require it.
	 *
	 *
	 * Iterating the vector in the way described above has the advantage that when we do computation on particles
	 * and its neighborhood with the sequence described above it will happen that:
	 *
	 * * If to process a particle A we read some neighborhood particles to process the next particle A+1
	 *   we will probably read most of the previous particles.
	 *
	 *
	 * In order to show in practice what happen we first show the graph when we do not reorder
	 *
	 * \htmlinclude Vector/4_reorder/no_reorder.html
	 *
	 * The measure has oscillation but we see an asymptotic behavior from 0.04 in the initial condition to
	 * 0.124 . Below we show what happen when we use iterator from the Cell list hilbert
	 *
	 * \htmlinclude Vector/4_reorder/comp_reord.html
	 *
	 * In cases where particles does not move or move very slowly consider to use data-reordering, because it can
	 * give **8-10% speedup**
	 *
	 * \see \ref e4_reo
	 *
	 *  ## Timers ##
	 *
	 *  In order to collect the time of the force calculation we insert two timers around the function
	 *  calc_force. The overall performance is instead calculated with another timer around the time stepping
	 *
	 * \snippet Vector/4_reorder/main_data_ord.cpp timer start
	 * \snippet Vector/4_reorder/main_data_ord.cpp timer stop
	 *
	 * \see \ref e3_md_vi
	 *
	 *
	 */

	//! \cond [md steps] \endcond

	// Get the Cell list structure
	auto NN = vd.getCellList_hilb(r_cut);

	// calculate forces
	calc_forces(vd,NN,sigma12,sigma6);
	unsigned long int f = 0;

	timer time2;
	time2.start();

	#ifndef TEST_RUN
	size_t Nstep = 30000;
	#else
	size_t Nstep = 300;
	#endif

	// MD time stepping
	for (size_t i = 0; i < Nstep ; i++)
	{
		// Get the iterator
		auto it3 = vd.getDomainIterator();

		// integrate velicity and space based on the calculated forces (Step1)
		while (it3.isNext())
		{
			auto p = it3.get();

			// here we calculate v(tn + 0.5)
			vd.template getProp<velocity>(p)[0] += 0.5*dt*vd.template getProp<force>(p)[0];
			vd.template getProp<velocity>(p)[1] += 0.5*dt*vd.template getProp<force>(p)[1];
			vd.template getProp<velocity>(p)[2] += 0.5*dt*vd.template getProp<force>(p)[2];

			// here we calculate x(tn + 1)
			vd.getPos(p)[0] += vd.template getProp<velocity>(p)[0]*dt;
			vd.getPos(p)[1] += vd.template getProp<velocity>(p)[1]*dt;
			vd.getPos(p)[2] += vd.template getProp<velocity>(p)[2]*dt;

			++it3;
		}

		// Because we mooved the particles in space we have to map them and re-sync the ghost
		vd.map();
		vd.template ghost_get<>();

		timer time;
		if (i % 10 == 0)
			time.start();

		// calculate forces or a(tn + 1) Step 2
		calc_forces(vd,NN,sigma12,sigma6);

		if (i % 10 == 0)
		{
			time.stop();
			x.add(i);
			y.add({time.getwct()});
		}

		// Integrate the velocity Step 3
		auto it4 = vd.getDomainIterator();

		while (it4.isNext())
		{
			auto p = it4.get();

			// here we calculate v(tn + 1)
			vd.template getProp<velocity>(p)[0] += 0.5*dt*vd.template getProp<force>(p)[0];
			vd.template getProp<velocity>(p)[1] += 0.5*dt*vd.template getProp<force>(p)[1];
			vd.template getProp<velocity>(p)[2] += 0.5*dt*vd.template getProp<force>(p)[2];

			++it4;
		}

		// After every iteration collect some statistic about the confoguration
		if (i % 100 == 0)
		{	
			// We write the particle position for visualization (Without ghost)
			vd.deleteGhost();
			vd.write("particles_",f);

			// we resync the ghost
			vd.ghost_get<>();

			// We calculate the energy
			double energy = calc_energy(vd,NN,sigma12,sigma6);
			auto & vcl = create_vcluster();
			vcl.sum(energy);
			vcl.execute();

			// We also print on terminal the value of the energy
			// only one processor (master) write on terminal
			if (vcl.getProcessUnitID() == 0)
				std::cout << std::endl << "Energy: " << energy << std::endl;

			f++;
		}
	}

	time2.stop();
	std::cout << "Performance: " << time2.getwct() << std::endl;

	//! \cond [md steps] \endcond

	/*!
	 * \page Vector_4_comp_reo Vector 4 computational reordering and cache friendliness
	 *
	 * ## Plotting graphs ##
	 *
	 * After we terminate the MD steps our vector x contains at which iteration we benchmark the force
	 * calculation time, while y contains the measured time at that time-step. We can produce a graph X Y
	 *
	 * \note The graph produced is an svg graph that can be view with a browser. From the browser we can
	 *       also easily save the graph into pure svg format
	 *
	 * \snippet Vector/4_reorder/main_comp_ord.cpp google chart
	 *
	 */

	//! \cond [google chart] \endcond

	// Google charts options, it store the options to draw the X Y graph
	GCoptions options;

	// Title of the graph
	options.title = std::string("Force calculation time");

	// Y axis name
	options.yAxis = std::string("Time");

	// X axis name
	options.xAxis = std::string("iteration");

	// width of the line
	options.lineWidth = 1.0;

	// Object that draw the X Y graph
	GoogleChart cg;

	// Add the graph
	// The graph that it produce is in svg format that can be opened on browser
	cg.AddLinesGraph(x,y,options);

	// Write into html format
	cg.write("gc_plot2_out.html");

	//! \cond [google chart] \endcond

	/*!
	 *
	 * \page Vector_4_comp_reo Vector 4 computational reordering and cache friendliness
	 *
	 * ## Finalize ##
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Vector/4_reorder/main_comp_ord.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 *
	 * \page Vector_4_comp_reo Vector 4 computational reordering and cache friendliness
	 *
	 * # Full code #
	 *
	 * \include Vector/4_reorder/main_comp_ord.cpp
	 *
	 */
}




