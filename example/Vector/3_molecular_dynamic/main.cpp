/*! \page Vector_3_md Vector 3 molecular dynamic
 *
 * \subpage Vector_3_md_dyn
 * \subpage Vector_3_md_vl
 *
 */

#include "Vector/vector_dist.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include "timer.hpp"

/*!
 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
 *
 * [TOC]
 *
 * # Molecular Dynamic with Lennard-Jones potential # {#e3_md}
 *
 * This example show a simple Lennard-Jones molecular dynamic simulation in a stable regime.
 * Particle feel each other by the potential.
 *
 * \f$ V(x_p,x_q) = 4( (\frac{\sigma}{r})^{12} - (\frac{\sigma}{r})^6  ) \f$
 *
 * ## Constants ##
 *
 * Here we define some useful constants
 *
 * \snippet Vector/3_molecular_dynamic/main.cpp constants
 *
 */

//! \cond [constants] \endcond

constexpr int velocity = 0;
constexpr int force = 1;

//! \cond [constants] \endcond

/*!
 *
 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
 *
 * ## Calculate forces ## {#e3_md_cf}
 *
 * In this function we calculate the forces between particles. It require the vector of particles
 * Cell list and scaling factor for the Lennard-Jhones potential.
 *
 * \snippet Vector/3_molecular_dynamic/main.cpp calc forces
 * \snippet Vector/3_molecular_dynamic/main.cpp calc forces2
 *
 *
 * In the following we are going into the detail of this function
 *
 */

//! \cond [calc forces] \endcond

template<typename CellList> void calc_forces(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, CellList & NN, double sigma12, double sigma6, double r_cut2)
{

//! \cond [calc forces] \endcond

	/*!
	 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
	 *
	 * This function in called several time and require neighborhood of each particle. In order to speed-up the
	 * Cell-list construction we can use updateCellList function to reuse the memory of the previous cell-list.
	 * updateCellList can be faster than createCellList
	 *
	 * \see \ref e1_part_celllist
	 *
	 * \snippet Vector/3_molecular_dynamic/main.cpp ucl
	 *
	 */

	//! \cond [ucl] \endcond

	vd.updateCellList(NN);

	//! \cond [ucl] \endcond

	/*!
	 *
	 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
	 *
	 * Get an iterator over the particles and get its position. For each particle p iterate in its neighborhood q
	 * and calculate the force based on the Lennard-Jhones potential given by
	 *
	 * \f$ F(x_p,x_q) = 24( \frac{2 \sigma^{12}}{r^{13}} - \frac{\sigma^6}{r^{7}}) \hat{r} \f$
	 *
	 * \see \ref e0_s_assign_pos
	 *
	 * \snippet Vector/3_molecular_dynamic/main.cpp force calc
	 *
	 */

	//! \cond [force calc] \endcond

	// Get an iterator over particles
	auto it2 = vd.getDomainIterator();

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
			if (q == p.getKey())	{++Np; continue;};

			// Get the position of p
			Point<3,double> xq = vd.getPos(q);

			// Get the distance between p and q
			Point<3,double> r = xp - xq;

			// take the norm of this vector
			double rn = norm2(r);

			if (rn > r_cut2)
			{++Np; continue;};

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

	//! \cond [force calc] \endcond

//! \cond [calc forces2] \endcond

}

//! \cond [calc forces2] \endcond

/*!
 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
 *
 * ## Calculate energy ## {#e3_md_ce}
 *
 * We also need a function to calculate energy, this function require the same parameter as calculate forces
 *
 * \snippet Vector/3_molecular_dynamic/main.cpp calc energy
 * \snippet Vector/3_molecular_dynamic/main.cpp calc energy2
 *
 */

//! \cond [calc energy] \endcond

template<typename CellList> double calc_energy(vector_dist<3,double, aggregate<double[3],double[3]> > & vd, CellList & NN, double sigma12, double sigma6, double r_cut2)
{
        double rc = r_cut2;
        double shift = 2.0 * ( sigma12 / (rc*rc*rc*rc*rc*rc) - sigma6 / ( rc*rc*rc) );

//! \cond [calc energy] \endcond

	/*!
	 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
	 *
	 * Reset the counter for the energy counter and
	 * update the cell list from the actual particle configuration
	 *
	 * \snippet Vector/3_molecular_dynamic/main.cpp up cell ene
	 *
	 * In the following we are going into the detail of this function
	 *
	 */

	//! \cond [up cell ene] \endcond

	double E = 0.0;
//	vd.updateCellList(NN);

	//! \cond [up cell ene] \endcond

	/*!
	 *
	 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
	 *
	 * First we get an iterator over the particles and get its position. For each particle p iterate in its neighborhood q
	 * and calculate the energy based on the Lennard-Jhones potential given by
	 *
	 * \f$ V(x_p,x_q) = 4(\frac{1}{r^{12}} - \frac{1}{r^{6}}) \f$
	 *
	 * \see \ref e0_s_assign_pos
	 *
	 * \snippet Vector/3_molecular_dynamic/main.cpp energy calc comp
	 *
	 */

	//! \cond [energy calc comp] \endcond

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
		auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(p)));

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

			if (rn > r_cut2)
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

	//! \cond [energy calc comp] \endcond

//! \cond [calc energy2] \endcond

	// Calculated energy
	return E;
}

//! \cond [calc energy2] \endcond

int main(int argc, char* argv[])
{
	/*!
	 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
	 *
	 * ## Initialization ## {#e3_md_init}
	 *
	 * After we defined the two main function calc forces and calc energy, we Initialize
	 *  the library, we create a Box that define our domain, boundary conditions and ghost.
	 *  Than we define important parameters of the simulation, time step integration,
	 * size of the box, and cut-off radius of the interaction. We also define 2 vectors
	 * x and y (they are like std::vector) used for statistic
	 *
	 * \see \ref e0_s_init
	 *
	 * \snippet Vector/3_molecular_dynamic/main.cpp constants run
	 *
	 */

	//! \cond [constants run] \endcond

	openfpm_init(&argc,&argv);

	double sigma = 0.1;
	double r_cut = 3.0*sigma;

	// we will use it do place particles on a 10x10x10 Grid like
	size_t sz[3] = {10,10,10};

	// domain
	Box<3,double> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,double> ghost(r_cut);

	double dt = 0.0005;
	double sigma12 = pow(sigma,12);
	double sigma6 = pow(sigma,6);

	openfpm::vector<double> x;
	openfpm::vector<openfpm::vector<double>> y;

	//! \cond [constants run] \endcond

	/*!
	 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
	 *
	 * Than we define a distributed vector in 3D, containing 2 vectorial properties the
	 * first is the actual velocity of the particle the other is the force
	 *
	 * \see \ref e0_s_vector_inst
	 *
	 * \snippet Vector/3_molecular_dynamic/main.cpp vect create
	 *
	 */

	//! \cond [vect create] \endcond

	vector_dist<3,double, aggregate<double[3],double[3]> > vd(0,box,bc,ghost);

	//! \cond [vect create] \endcond

	/*!
	 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
	 *
	 * ## Particles on a grid like position ## {#e3_md_gl}
	 *
	 * We define a grid iterator, to create particles on a grid like way. In the same cycle we also reset
	 * force and velocity
	 *
	 * \see \ref e1_cl_gr_it
	 *
	 * \snippet Vector/3_molecular_dynamic/main.cpp vect grid
	 *
	 */

	//! \cond [vect grid] \endcond

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

	vd.map();
	vd.ghost_get<>();

	//! \cond [vect grid] \endcond

	/*!
	 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
	 *
	 * ## Molecular dynamic steps ## {#e3_md_vi}
	 *
	 * Here we do 10000 MD steps using verlet integrator
	 *
	 * The verlet integration stepping look like this
	 *
	 * \f[ \vec{v}(t_{n}+1/2) = \vec{v}_p(t_n) + \frac{1}{2} \delta t \vec{a}(t_n) \f]
	 * \f[ \vec{x}(t_{n}+1) = \vec{x}_p(t_n) + \delta t \vec{v}(t_n+1/2) \f]
	 *
	 * calculate the forces from \f$ \vec{a} (t_{n}) \f$ finally
	 *
	 * \f[ \vec{v}(t_{n+1}) = \vec{v}_p(t_n+1/2) + \frac{1}{2} \delta t \vec{a}(t_n+1) \f]
	 *
	 * The cell-list structure is required to calculate forces. As demonstration
	 * purpose instead of using the standard Cell-list with (getCellList) we use the CELL_MEMBAL
	 * type. The impact in performance of using the CELL_MEMBAL instead of CELL_MEMFAST is less
	 * than 1% on the other hand CELL_MEMFAST use 16 Megabyte of memory
	 *
	 * \see \ref e1_cls_types
	 *
	 * Inside this cycle we are using several features that has been explained before in particular
	 *
	 * \see \ref e0_s_assign_pos
	 *
	 * \see \ref e0_s_map
	 *
	 * \see \ref e0_s_reduce
	 *
	 * \see \ref e1_part_ghost
	 *
	 * \see \ref e0_s_vis_vtk
	 *
	 * \snippet Vector/3_molecular_dynamic/main.cpp md steps
	 *
	 */

	timer tsim;
	tsim.start();

	//! \cond [md steps] \endcond

	// Get the Cell list structure
	auto NN = vd.getCellList<CELL_MEMBAL<3,double>>(r_cut);

	// The standard
	// auto NN = vd.getCellList(r_cut);

	// calculate forces
	calc_forces(vd,NN,sigma12,sigma6,r_cut*r_cut);
	unsigned long int f = 0;

	// MD time stepping
	for (size_t i = 0; i < 10000 ; i++)
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

		// Because we moved the particles in space we have to map them and re-sync the ghost
		vd.map();
		vd.template ghost_get<>();

		// calculate forces or a(tn + 1) Step 2
		calc_forces(vd,NN,sigma12,sigma6,r_cut*r_cut);

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

                // After every iteration collect some statistic about the configuration
                if (i % 100 == 0)
                {
                        // We write the particle position for visualization (Without ghost)
                        vd.deleteGhost();
                        vd.write_frame("particles_",f);

                        // we resync the ghost
                        vd.ghost_get<>();

                        // We calculate the energy
                        double energy = calc_energy(vd,NN,sigma12,sigma6,r_cut*r_cut);
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

	/*!
	 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
	 *
	 * ## Plotting graphs ## {#e3_md_pg}
	 *
	 * After we terminate the MD steps our vector x contains at which iteration we calculated the energy
	 * while y contains the energy value at that time-step. We can produce a graph X Y
	 *
	 * \note The graph produced is an svg graph that can be view with a browser. From the browser we can
	 *       also easily save the graph into pure svg format
	 *
	 * \snippet Vector/3_molecular_dynamic/main.cpp google chart
	 *
	 */

	//! \cond [google chart] \endcond

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

	//! \cond [google chart] \endcond

	/*!
	 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
	 *
	 * ## Finalize ## {#finalize_v_e3_md}
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Vector/3_molecular_dynamic/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
	 *
	 * ## Full code ## {#code_v_e3_md}
	 *
	 * \include Vector/3_molecular_dynamic/main.cpp
	 *
	 */

	/*!
	 * \page Vector_3_md_dyn Vector 3 molecular dynamic with cell-list
	 *
	 * ## Code with expression ## {#code_v_e3_md_expr}
	 *
	 * Here we also show how we can simplify the example using expressions
	 *
	 * \see \ref Vector_2_expression
	 *
	 * \include Vector/3_molecular_dynamic/main_expr.cpp
	 *
	 */
}




