
#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include "timer.hpp"

/*!
 * \page Vector_3_md_vl Vector 3 molecular dynamic with Verlet list
 *
 *
 * [TOC]
 *
 * # Molecular Dynamic with Lennard-Jones potential with verlet list # {#e3_md_vl}
 *
 * This example show a simple Lennard-Jones molecular dynamic simulation in a stable regime.
 * We will use Verlet-list in order to get a speed-up from force calculation
 *
 * ## Constants ##
 *
 * Here we define some useful constants
 *
 * \snippet Vector/3_molecular_dynamic/main_vl.cpp constants
 *
 */

//! \cond [constants] \endcond

constexpr int velocity = 0;
constexpr int force = 1;

//! \cond [constants] \endcond

/*!
 *
 * \page Vector_3_md_vl Vector 3 molecular dynamic with Verlet list
 *
 * ## Calculate forces ## {#e3_md_vl_cf}
 *
 * In this function we calculate the forces between particles. It require the vector of particles,
 * the Verlet-list and sigma for the Lennard-Jhones potential. The function is exactly the same
 * as the original with the following changes
 *
 * \see \ref e3_md_cf
 *
 * * The function accept a VerletList instead of a CellList
 *  \snippet main_vl.cpp arg diff
 *
 * * There is no call to updateCellList
 *
 * * How to get an iterator over neighborhood of a particle
 *   \snippet main_vl.cpp NN iterator
 *
 * Teh rest remain the same
 *
 * \snippet Vector/3_molecular_dynamic/main_vl.cpp calc forces vl
 *
 */

//! \cond [calc forces vl] \endcond

//! \cond [arg diff] \endcond

void calc_forces(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, VerletList<3, double, VL_NON_SYMMETRIC, Mem_fast<>, shift<3, double> > & NN, double sigma12, double sigma6, double r_cut)
{
	//! \cond [arg diff] \endcond

	// Get an iterator over particles
	auto it2 = vd.getDomainIterator();

	// For each particle p ...
	while (it2.isNext())
	{
		// ... get the particle p
		auto p = it2.get();

		// Get the position xp of the particle
		Point<3,double> xp = vd.getPos(p);

		// Reset the forice counter
		vd.template getProp<force>(p)[0] = 0.0;
		vd.template getProp<force>(p)[1] = 0.0;
		vd.template getProp<force>(p)[2] = 0.0;

		//! \cond [NN iterator] \endcond

		// Get an iterator over the neighborhood particles of p
		auto Np = NN.getNNIterator(p.getKey());

		//! \cond [NN iterator] \endcond

		// For each neighborhood particle ...
		while (Np.isNext())
		{
			// ... q
			auto q = Np.get();

			// if (p == q) skip this particle
			if (q == p.getKey())	{++Np; continue;};

			// Get the position of p
			Point<3,double> xq = vd.getPos(q);

			// Get the distance between p and q
			Point<3,double> r = xp - xq;

			// take the norm of this vector
			double rn = norm2(r);

			if (rn > r_cut * r_cut) {++Np;continue;}

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

//! \cond [calc forces vl] \endcond

/*!
 * \page Vector_3_md_vl Vector 3 molecular dynamic with verlet list
 *
 * ## Calculate energy ## {#e3_md_vl_ce}
 *
 * We also need a function to calculate energy, this function require the same parameter as calculate forces
 *
 * \see \ref e3_md_ce
 *
 * The following changes has been made
 *
 * * The function accept a VerletList instead of a cell-List
 * * There is no call to updateCellList
 * * How to get an iterator over neigborhood particles
 *
 * \snippet Vector/3_molecular_dynamic/main_vl.cpp calc energy vl
 *
 */

//! \cond [calc energy vl] \endcond

double calc_energy(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, VerletList<3, double, VL_NON_SYMMETRIC, Mem_fast<>, shift<3, double> > & NN, double sigma12, double sigma6, double r_cut)
{
	double E = 0.0;

	double rc = r_cut*r_cut;
	double shift = 2.0 * ( sigma12 / (rc*rc*rc*rc*rc*rc) - sigma6 / ( rc*rc*rc) );

	// Get the iterator
	auto it2 = vd.getDomainIterator();

	// For each particle ...
	while (it2.isNext())
	{
		// ... p
		auto p = it2.get();

		// Get the position of the particle p
		Point<3,double> xp = vd.getPos(p);

		// Get an iterator over the neighborhood of the particle p
		auto Np = NN.getNNIterator(p.getKey());

		// For each neighborhood of the particle p
		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			// if p == q skip this particle
			if (q == p.getKey())	{++Np; continue;};

			// Get position of the particle q
			Point<3,double> xq = vd.getPos(q);

			// take the normalized direction
			double rn = norm2(xp - xq);

			if (rn >= r_cut*r_cut)
			{++Np;continue;}

			// potential energy (using pow is slower)
			E += 2.0 * ( sigma12 / (rn*rn*rn*rn*rn*rn) - sigma6 / ( rn*rn*rn) ) - shift;

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

//! \cond [calc energy vl] \endcond

int main(int argc, char* argv[])
{
	/*!
	 * \page Vector_3_md_vl Vector 3 molecular dynamic with Verlet list
	 *
	 * ## Simulation ## {#e3_md_vl_sim}
	 *
	 * The simulation is equal to the simulation explained in the example molecular dynamic
	 *
	 * \see \ref e3_md
	 *
	 * The differences are that:
	 *
	 * * The Ghost must be r_cut+skin
	 * \snippet
	 *
	 * * We create a Verlet list with skin instead of a Cell list
	 * \snippet Vector/3_molecular_dynamic/main_vl.cpp verlet skin
	 *
	 * * every 10 steps we do a map and update the verlet-list, in all the other case we just do skip labelling
	 * \snippet Vector/3_molecular_dynamic/main_vl.cpp update verlet
	 *
	 * **Explanation**
	 *
	 * Updating the verlet list is extremely expensive. For this reason we create a Verlet list
	 * that contain r_cut + skin particles. Using the fact that during the full simulation each
	 *  particle does not move more than 0.0015 in one iteration, if the skin is 0.03
	 *  we can update the Verlet list every \f$ \frac{0.03}{2 \cdot 0.0015} = 10 \f$. The 2 factor if given by the fact
	 *  that in the worst case where one particle is going left and one on the right from the prospective of
	 *  one particle the particle moove \f$ 2 \cdot 0.0015 \f$.
	 *
	 * Because the Verlet lists are constructed based on the local-id of the particles a map or a ghost_get
	 *  would invalidate the verlet. For this reason the map is called every 10 time-step (when we
	 *  update the verlet), and a particular ghost_get with SKIP_LABELLING is used during every iteration.
	 *
	 *  The function ghost_get with skip labeling does not recompute the particle to send but use the
	 *  the ids of the old particles updating the positions (and properties if needed) and keeping the old
	 *   indexes without invalidating the Verlet-list. Doing this we can avoid to send particles that are
	 *   entering the ghost area r_cut+skin. Because we know that no particle in 10 iteration can travel for a
	 *   distance bigger than the skin, we are sure that in 10 iteration no-new particle that were not in the
	 *   r_cut+skin ghost area can enter the ghost area r_cut.
	 *
	 *
	 * \snippet Vector/3_molecular_dynamic/main_vl.cpp simulation
	 *
	 */

	//! \cond [simulation] \endcond

	double dt = 0.00025;
	double sigma = 0.1;
	double r_cut = 3.0*sigma;
	double r_gskin = 1.3*r_cut;
	double sigma12 = pow(sigma,12);
	double sigma6 = pow(sigma,6);

	openfpm::vector<double> x;
	openfpm::vector<openfpm::vector<double>> y;

	openfpm_init(&argc,&argv);
	Vcluster<> & v_cl = create_vcluster();

	// we will use it do place particles on a 10x10x10 Grid like
	size_t sz[3] = {10,10,10};

	// domain
	Box<3,double> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,double> ghost(r_gskin);

	vector_dist<3,double, aggregate<double[3],double[3]> > vd(0,box,bc,ghost);

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

	timer tsim;
	tsim.start();

	//! \cond [verlet skin] \endcond

	// Get the Cell list structure
	auto NN = vd.getVerlet(r_gskin);

	//! \cond [verlet skin] \endcond

	// calculate forces
	calc_forces(vd,NN,sigma12,sigma6,r_cut);
	unsigned long int f = 0;

	int cnt = 0;
	double max_disp = 0.0;

	// MD time stepping
	for (size_t i = 0; i < 10000 ; i++)
	{
		// Get the iterator
		auto it3 = vd.getDomainIterator();

		double max_displ = 0.0;

		// integrate velicity and space based on the calculated forces (Step1)
		while (it3.isNext())
		{
			auto p = it3.get();

			// here we calculate v(tn + 0.5)
			vd.template getProp<velocity>(p)[0] += 0.5*dt*vd.template getProp<force>(p)[0];
			vd.template getProp<velocity>(p)[1] += 0.5*dt*vd.template getProp<force>(p)[1];
			vd.template getProp<velocity>(p)[2] += 0.5*dt*vd.template getProp<force>(p)[2];

			Point<3,double> disp({vd.template getProp<velocity>(p)[0]*dt,vd.template getProp<velocity>(p)[1]*dt,vd.template getProp<velocity>(p)[2]*dt});

			// here we calculate x(tn + 1)
			vd.getPos(p)[0] += disp.get(0);
			vd.getPos(p)[1] += disp.get(1);
			vd.getPos(p)[2] += disp.get(2);

			if (disp.norm() > max_displ)
				max_displ = disp.norm();

			++it3;
		}

		if (max_disp < max_displ)
			max_disp = max_displ;

		//! \cond [update verlet] \endcond

		// Because we moved the particles in space we have to map them and re-sync the ghost
		if (cnt % 10 == 0)
		{
			vd.map();
			vd.template ghost_get<>();
			// Get the Cell list structure
			vd.updateVerlet(NN,r_gskin);
		}
		else
		{
			vd.template ghost_get<>(SKIP_LABELLING);
		}

		//! \cond [update verlet] \endcond

		cnt++;

		// calculate forces or a(tn + 1) Step 2
		calc_forces(vd,NN,sigma12,sigma6,r_cut);

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
			vd.write_frame("particles_",f);

			// we resync the ghost
			vd.ghost_get<>(SKIP_LABELLING);

			// We calculate the energy
			double energy = calc_energy(vd,NN,sigma12,sigma6,r_cut);
			auto & vcl = create_vcluster();
			vcl.sum(energy);
			vcl.max(max_disp);
			vcl.execute();

			// we save the energy calculated at time step i c contain the time-step y contain the energy
			x.add(i);
			y.add({energy});

			// We also print on terminal the value of the energy
			// only one processor (master) write on terminal
			if (vcl.getProcessUnitID() == 0)
				std::cout << "Energy: " << energy << "   " << max_disp << "  " << std::endl;

			max_disp = 0.0;

			f++;
		}
	}

	tsim.stop();
	std::cout << "Time: " << tsim.getwct()  << std::endl;

	//! \cond [simulation] \endcond

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

	// Object that draw the X Y graph
	GoogleChart cg;

	// Add the graph
	// The graph that it produce is in svg format that can be opened on browser
	cg.AddLinesGraph(x,y,options);

	// Write into html format
	cg.write("gc_plot2_out.html");

	//! \cond [google chart] \endcond

	/*!
	 * \page Vector_3_md_vl Vector 3 molecular dynamic with Verlet list
	 *
	 * ## Finalize ## {#finalize_v_e3_md_vl}
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Vector/3_molecular_dynamic/main_vl.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond
}
