#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include "timer.hpp"

/*!
 * \page Vector_5_md_vl_sym Vector 5 molecular dynamic with symmetric Verlet list
 *
 *
 * # Molecular dynamic with symmetric interactions # {#md_e5_sym}
 *
 * In a previous example we show how to build a parallel molecular dynamic simulation.
 * We also show how it was possible to achive better performance using Verlet list.
 *
 * \see \ref e3_md_vl
 *
 *  In this example we show how to improve even more our code using symmetric Verlet-list.
 *  If we look at the form of the Lennard-Jhones potential we will see that
 * calculating the force that the particle **p** produce on **q** is equivalent to the force
 * that **q** produce on **p**. This mean that when we use normal verlet-list we are redundantly
 * doing double calculation. In order to avoid it we will use symmetric verlet-list. Symmetric
 * verlet list have their feature in store the pair **p,q** only one time. If **p** store **q**
 * as neighborhood than **q** does not store **p** as neighborhood. This Mean that when we calculate
 * the contribution for of **q** to **p** we have to add such contribution also to **q**.
 *
 * The example is exactly equivalent to the non-symmetric with few differences in **calc_forces** and
 * **calc_energies**
 *
 *
 */

//! \cond [constants] \endcond

constexpr int velocity = 0;
constexpr int force = 1;

//! \cond [constants] \endcond

/*!
 *
 * \page Vector_5_md_vl_sym Vector 5 molecular dynamic with symmetric Verlet list
 *
 * ## Calculate forces ## {#md_e5_calc_force}
 *
 * In this function we calculate the forces between particles. It require the vector of particles
 * Cell list and sigma factor for the Lennard-Jhones potential. The function is exactly the same
 * as the original
 *
 * \see \ref e3_md_vl_cf
 *
 * with the following changes
 *
 *
 * * If we calculate the force for **p-q** we are also adding this force to **q-p**
 *  \snippet Vector/5_molecular_dynamic_sym/main.cpp add to q
 *
 * * At the end of the calculation we have to execute a ghost put
 *   \snippet Vector/5_molecular_dynamic_sym/main.cpp ghost_put
 *
 * ###Explanation### {#md_e5_sym_expl}
 *
 * The first point is given by the fact that if the pair is stored once, when we calculate
 * the force, we have to add the contribution to both particles. The second instead is
 * is given by the fact that **q** can be a ghost particles. In case **q** is a ghost particle
 * we are adding the contribution to the ghost particle and not to the real one. To
 * add the contribution to the real particle we have to use the function **ghost_put**.
 * This function send back the information to the original processor, that will merge the
 * information (in this case add_)
 *
 *
 * \snippet Vector/5_molecular_dynamic_sym/main.cpp calc forces vl
 *
 */

//! \cond [calc forces vl] \endcond

template<typename VerletList>
void calc_forces(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, VerletList & NN, double sigma12, double sigma6, double r_cut)
{
	// Reset force on the ghost

	auto itg = vd.getDomainAndGhostIterator();

	while (itg.isNext())
	{
		auto p = itg.get();

		// Reset force
		vd.getProp<force>(p)[0] = 0.0;
		vd.getProp<force>(p)[1] = 0.0;
		vd.getProp<force>(p)[2] = 0.0;

		++itg;
	}

	//! \cond [real and ghost] \endcond

	// Get an iterator over particles
	auto it2 = vd.getDomainIterator();

	//! \cond [real and ghost] \endcond

	// For each particle p ...
	while (it2.isNext())
	{
		// ... get the particle p
		auto p = it2.get();

		// Get the position xp of the particle
		Point<3,double> xp = vd.getPos(p);

		// Get an iterator over the neighborhood particles of p
		// Note that in case of symmetric
		auto Np = NN.getNNIterator(p.getKey());

		// For each neighborhood particle ...
		while (Np.isNext())
		{
			// ... q
			auto q = Np.get();

			// if (p == q) skip this particle
			if (q == p.getKey())	{++Np; continue;};

			// Get the position of q
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

			//! \cond [add to q] \endcond

			// we sum the force produced by p on q
			vd.template getProp<force>(q)[0] -= f.get(0);
			vd.template getProp<force>(q)[1] -= f.get(1);
			vd.template getProp<force>(q)[2] -= f.get(2);

			//! \cond [add to q] \endcond

			// Next neighborhood
			++Np;
		}

		// Next particle
		++it2;
	}

	//! \cond [ghost_put] \endcond

	// Sum the contribution to the real particles
	vd.ghost_put<add_,force>();

	//! \cond [ghost_put] \endcond
}

//! \cond [calc forces vl] \endcond


/*!
 *
 * \page Vector_5_md_vl_sym Vector 5 molecular dynamic with symmetric Verlet list
 *
 * ## Calculate energy ## {#md_e5_calc_ene}
 *
 * For the energy we use symmetric verlet-list in the same way as we did for calc_forces.
 * Because the symmetric verlet-list span each couple one time, we have to remove the division by two
 * (in this case we use the original factor 4.0 of the Lennard-Jhones potential rather than 2.0).
 *
 * \snippet Vector/5_molecular_dynamic_sym/main.cpp calc energy vl
 *
 */

//! \cond [calc energy vl] \endcond

template<typename VerletList>
double calc_energy(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, VerletList & NN, double sigma12, double sigma6, double r_cut)
{
	double E = 0.0;

	double rc = r_cut*r_cut;
	double shift = 4.0 * ( sigma12 / (rc*rc*rc*rc*rc*rc) - sigma6 / ( rc*rc*rc) );

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

		double Ep = E;

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

			if (rn >= r_cut*r_cut)	{++Np;continue;}

			// potential energy (using pow is slower)
			E += 4.0 * ( sigma12 / (rn*rn*rn*rn*rn*rn) - sigma6 / ( rn*rn*rn) ) - shift;

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
	 * \page Vector_5_md_vl_sym Vector 5 molecular dynamic with symmetric Verlet list
	 *
	 * ## Simulation ## {#md_e5_sym_sim}
	 *
	 * The simulation is equal to the simulation explained in the example molecular dynamic
	 *
	 * \see \ref e3_md_vl
	 *
	 * The difference is that we create a symmetric Verlet-list instead of a normal one
	 * \snippet Vector/5_molecular_dynamic_sym/main.cpp sim verlet
	 *
	 * The rest of the code remain unchanged
	 *
	 * \snippet Vector/5_molecular_dynamic_sym/main.cpp simulation
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

	size_t k = 0;
	size_t start = vd.accum();

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

		k++;
		++it;
	}

	vd.map();
	vd.ghost_get<>();

	timer tsim;
	tsim.start();

	//! \cond [sim verlet] \endcond

	// Get the Cell list structure
	auto NN = vd.getVerletSym(r_gskin);

	//! \cond [sim verlet] \endcond

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
			vd.write("particles_",f);

			// we resync the ghost
			vd.ghost_get<>();


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
				std::cout << "Energy: " << energy << "   " << max_disp << "   " << std::endl;

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
	 * \page Vector_5_md_vl_sym Vector 5 molecular dynamic with symmetric Verlet list
	 *
	 * ## Finalize ## {#finalize_v_e5_md_sym}
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Vector/5_molecular_dynamic_sym/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_5_md_vl_sym Vector 5 molecular dynamic with symmetric Verlet list
	 *
	 * ## Full code ## {#full_code_v_e5_md_sym}
	 *
	 * \include Vector/5_molecular_dynamic_sym/main.cpp
	 *
	 */
}
