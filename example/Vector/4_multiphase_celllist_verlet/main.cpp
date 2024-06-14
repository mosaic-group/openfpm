
#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "NN/CellList/multiphase/CellListM.hpp"
#include "Vector/vector_dist_multiphase_functions.hpp"

/*!
 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list and verlet
 *
 * [TOC]
 *
 *
 * # Vector Multi Phase cell-list and Verlet # {#e4_ph_cl}
 *
 * This example show how to use multi-phase cell lists and Verlet-list.More in general
 * it show how to construct Verlet and Cell-list between multiple vector_dist.
 *
 *
 */

int main(int argc, char* argv[])
{
	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * ## Initialization ##
	 *
	 * Here we Initialize the library, and we create a set of distributed vectors all forced to have the same
	 * decomposition. Each vector identify one phase.
	 *
	 * \note Be carefull on how you initialize the other phases. All the other phases
	 *       must be forced to use the same decomposition. In order to do this we have
	 *       to use the special constructor where we pass the decomposition from the first
	 *       phase. The second parameter is just the the number of particles
	 *
	 * \snippet Vector/4_multiphase_celllist_verlet/main.cpp Initialization and parameters
	 *
	 */

	//! \cond [Initialization and parameters] \endcond

	openfpm_init(&argc,&argv);

	// Vcluster for general usage
	Vcluster<> & v_cl = create_vcluster();

	// The domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// cut-off radius
	float r_cut = 0.05;

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(r_cut);

	// The set of phases
	openfpm::vector< vector_dist<3,float, aggregate<double,double>> > phases;

	// first phase
	phases.add( vector_dist<3,float, aggregate<double,double>>(4096,box,bc,ghost) );

	// The other 3 phases
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(0).getDecomposition(),4096) );
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(0).getDecomposition(),4096) );
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(0).getDecomposition(),4096) );


	//! \cond [Initialization and parameters] \endcond


	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * ## Create phases ##
	 *
	 * We initialize all the phases with particle randomly positioned in the space. Completed this
	 * iteration we redistribute the particles using the classical map, and we synchronize the ghost
	 *
	 * \warning we cannot use the same iterator for all the phases even if the number of particles for
	 *          each phase is the same. Consider that the number of particles (per-processor) can be
	 *          different. The case of the initialization ("before map") is the only case where we are
	 *          sure that the number of particles per processor is the same for each phase
	 *
	 * \snippet Vector/4_multiphase_celllist_verlet/main.cpp rand dist
	 *
	 */

	//! \cond [rand dist] \endcond

	// An iterator over the particles of the phase0
	auto it = phases.get(0).getDomainIterator();

	// For all the particles of phase0
	while (it.isNext())
	{
		// particle p
		auto p = it.get();

		// Assign the position of the particles to each phase
		for (size_t i = 0 ; i < phases.size(); i++)
		{
			// Assign random position
			phases.get(i).getPos(p)[0] = (float)rand() / RAND_MAX;
			phases.get(i).getPos(p)[1] = (float)rand() / RAND_MAX;
			phases.get(i).getPos(p)[2] = (float)rand() / RAND_MAX;
		}

		// Next particle
		++it;
	}


	// Redistribute and sync all the phases
	for (size_t i = 0 ; i < 4 ; i++)
	{
		phases.get(i).map();
		phases.get(i).ghost_get<>();
	}

	//! \cond [rand dist] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * ## Multi phase Verlet-list ##
	 *
	 * In general verlet-list can be constructed from the vector itself using **getVerlet()**
	 *
	 * \see \ref Vector_3_md_vl
	 *
	 * In the multi-phase case if we use such function on phase0 such function
	 * produce a Verlet list where for each particle of the phase0 we get the neighborhood
	 *  particles within the phase0. Most of time we also want to construct Verlet-list across
	 *  phases, for example between phase0 and phase1. Let suppose now that we want to construct
	 *  a verlet-list that given the particles of the phase0 it return the neighborhood particles
	 *   in phase1. In order to do this we can create a Cell-list from phase1. Once we have the
	 *   Cell-list for phase 1 we give to **createVerlet()**
	 *  the particle set of phase0, the Cell-list of phase1 and the cut-off radius.
	 *  The function return a VerletList between phase0 and 1 that can be used to do computation.
	 *
	 *
	 * \snippet Vector/4_multiphase_celllist_verlet/main.cpp create multi-phase verlet
	 *
	 */

	//! \cond [create multi-phase verlet] \endcond

	{

	// Get the cell list of the phase1
	auto CL_phase1 = phases.get(1).getCellList(r_cut);

	// This function create a Verlet-list between phases 0 and 1
	auto NN_ver01 = createVerlet(phases.get(0),phases.get(1),CL_phase1,r_cut);

	//! \cond [create multi-phase verlet] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * Once we have a Verlet-list, we can do easily computation with that. We can use
	 * the function **getNNIterator()** to get an iterator of the neighborhood particles
	 *  for a specified particle. In this case for each particle of the phase0 we count
	 *  the neighborhood particles in phase1
	 *
	 * \snippet Vector/4_multiphase_celllist_verlet/main.cpp count part from phase0 to 1
	 *
	 */

	//! \cond [count part from phase0 to 1] \endcond

	// Get an iterator of the particles of the phase0
	it = phases.get(0).getDomainIterator();

	// For each particle of the phase0
	while (it.isNext())
	{
		// Get the particle p
		auto p = it.get();

		// reset the counter
		phases.get(0).getProp<0>(p) = 0;

		// Get an iterator over the neighborhood particles for the particle p
		auto Np = NN_ver01.getNNIterator(p.getKey());

		// For each neighborhood of the particle p
		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			// Count the number of particles
			// Here in general we can do our computation
			phases.get(0).getProp<0>(p)++;

			// Next particle
			++Np;
		}

		// Next particle
		++it;
	}

	}

	//! \cond [count part from phase0 to 1] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * ## From one phase to all ##
	 *
	 * It is also possible to construct a Verlet-list from one phase to all the other phases.
	 * To do this we use the function **createCellListM<2>()** to first create a
	 * multi-phase cell-list. The template parameter is required to indicate how many bit
	 * to reserve for the phase information.
	 *
	 * \note In case of
	 * 	* 2 bit mean that you can store up to 4 phases (2^62 number of particles for each phase)
	 * 	* 3 bit mean that you can store up to 8 phases (2^61 number of particles for each phase)
	 *
	 * Once created the multiphase-cell list from the phases. We can use the function
	 * **createVerletM()** to create a verlet-list from one phase to all the others.
	 * The function requires the particle for which we are constructing the verlet list, in this case
	 * the phase0, the Cell-list containing all the other phases and the cut-off radius.
	 *
	 *
	 * \snippet Vector/4_multiphase_celllist_verlet/main.cpp create multi-phase multi verlet
	 *
	 */

	//! \cond [create multi-phase multi verlet] \endcond

	// This function create an "Empty" Multiphase Cell List from all the phases

	// We just create a cell list with slightly bigget r_cut to avoid error in Verlet creation
	// because of round off errors
	float r_cut2 = r_cut*1.00001;
	auto CL_all = createCellListM<2>(phases,r_cut2);

	// This create a Verlet-list between phase0 and all the other phases
	auto NNver0_all = createVerletM<2>(0,phases.get(0),phases,CL_all,r_cut);

	//! \cond [create multi-phase multi verlet] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * Compute on a multiphase verlet-list is very similar to a non multi-phase.
	 * The only difference is that the function **get()** is substituted by two
	 * other functions **getP()** and **getV()** The first return the phase from where
	 * it come the particle the second it return the particle-id
	 *
	 * \snippet Vector/4_multiphase_celllist_verlet/main.cpp compute multi-phase multi verlet
	 *
	 */

	//! \cond [compute multi-phase multi verlet] \endcond

	// Get an iterator for the phase0 particles
	it = phases.get(0).getDomainIterator();

	while (it.isNext())
	{
		// Get the particle p
		auto p = it.get();

		// Get an interator over the neighborhood of the particle p
		auto Np = NNver0_all.getNNIterator(p.getKey());

		// reset the counter
		phases.get(0).getProp<0>(p) = 0;

		// For each neighborhood of the particle p
		while (Np.isNext())
		{
			// Get the particle q near to p
			auto q = Np.getP();

			// Get from which phase it come from
			auto ph_q = Np.getV();

			// count
			phases.get(0).getProp<0>(p)++;

			// Next particle
			++Np;
		}

		// Next particle
		++it;
	}

	//! \cond [compute multi-phase multi verlet] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * ## Symmetric interaction case ##
	 *
	 * The same functions exist also in the case we want to construct
	 * symmetric-verlet list.
	 *
	 * \see \ref Vector_5_md_vl_sym For more details on how to use symmetric
	 * verlet-list in a real case.
	 *
	 * In general the main differences can be summarized in
	 * * we have to reset the **forces** or **interaction** variable(s)
	 * * calculate the interaction p-q and apply the result to p and q at the same time
	 * * merge with ghost_put the result is stored in the ghost part with the
	 *   real particles
	 *
	 * \snippet Vector/4_multiphase_celllist_verlet/main.cpp compute sym multi-phase two phase
	 *
	 */

	{
	//! \cond [compute sym multi-phase two phase] \endcond

	// Get the cell list of the phase1
	auto CL_phase1 = phases.get(1).getCellListSym(r_cut);

	// This function create a Verlet-list between phases 0 and 1
	auto NN_ver01 = createVerletSym(phases.get(0),phases.get(1),CL_phase1,r_cut);

	// Get an iterator over the real and ghost particles
	it = phases.get(0).getDomainAndGhostIterator();

	// For each particles
	while (it.isNext())
	{
		// Get the particle p
		auto p = it.get();

		// Reset the counter
		phases.get(0).getProp<0>(p) = 0;

		// Next particle
		++it;
	}

	// Compute interaction from phase0 to phase1

	// Get an iterator over the real particles of phase0
	it = phases.get(0).getDomainIterator();

	// For each particle of the phase0
	while (it.isNext())
	{
		// get particle p
		auto p = it.get();

		// Get the neighborhood of particle p
		auto Np = NN_ver01.getNNIterator(p.getKey());

		// For each neighborhood of the particle p
		while (Np.isNext())
		{
			// Neighborhood particle q of p
			auto q = Np.get();

			// increment the counter for the phase 0 and 1
			phases.get(0).getProp<0>(p)++;
			phases.get(1).getProp<0>(q)++;

			// Next particle
			++Np;
		}

		// Next particle
		++it;
	}

	// Merge the information of the ghost with the real particles
	phases.get(0).ghost_put<add_,0>();
	phases.get(1).ghost_put<add_,0>();

	}

	//! \cond [compute sym multi-phase two phase] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * For the case of a Verlet-list between the phase0 and all the other phases.
	 * The code remain very similar to the previous one the only differences is that
	 * we have to create a Multi-phase cell list from the vector of the phases.
	 * When we compute with the verlet-list list now the simple function **get**
	 * is substituted by the function **getP** and **getV**. In this example we are
	 * creating 4 multiphase symmetric verlet-list
	 *
	 * * 0 to all
	 * * 1 to all
	 * * 2 to all
	 * * 3 to all
	 *
	 * The computation is an all to all
	 *
	 * \snippet Vector/4_multiphase_celllist_verlet/main.cpp create sym multi-phase multi verlet
	 *
	 */

	//! \cond [create sym multi-phase multi verlet] \endcond

	// Get an iterator over the phase0
	it = phases.get(0).getDomainAndGhostIterator();

	// For each particle of the phase 0
	while (it.isNext())
	{
		// get the particle p
		auto p = it.get();

		// Reset the counter for the domain and ghost particles
		phases.get(0).getProp<0>(p) = 0;
		phases.get(1).getProp<0>(p) = 0;
		phases.get(2).getProp<0>(p) = 0;
		phases.get(3).getProp<0>(p) = 0;

		// next particle
		++it;
	}

	// This function create an "Empty" Multiphase Cell List
	CL_all = createCellListSymM<2>(phases,r_cut);

	// Type of the multiphase Verlet-list
	typedef decltype(createVerletSymM<2>(0,phases.get(0),phases,CL_all,r_cut)) verlet_type;

	// for each phase we create one Verlet-list that contain the neighborhood
	// from all the phases
	verlet_type NNver_all[4];

	// Here we create a Verlet-list between each phases

	// 0 to all
	NNver_all[0] = createVerletSymM<2>(0,phases.get(0),phases,CL_all,r_cut);

	// 1 to all
	NNver_all[1] = createVerletSymM<2>(1,phases.get(1),phases,CL_all,r_cut);

	// 2 to all
	NNver_all[2] = createVerletSymM<2>(2,phases.get(2),phases,CL_all,r_cut);

	// 3 to all
	NNver_all[3] = createVerletSymM<2>(3,phases.get(3),phases,CL_all,r_cut);

	// For each phase
	for (size_t i = 0 ; i < phases.size() ; i++)
	{
		// Get an iterator over the particles of that phase
		it = phases.get(i).getDomainIterator();

		// for each particle of the phase
		while (it.isNext())
		{
			// Get the particle p
			auto p = it.get();

			// Get an iterator for neighborhood of the particle p
			auto Np = NNver_all[i].getNNIterator(p.getKey());

			// For each neighborhood particle
			while (Np.isNext())
			{
				// Get the particle q near to p
				auto q = Np.getP();

				// Get from which phase it come from
				auto ph_q = Np.getV();

				// increment the counter on both p and q
				phases.get(i).getProp<0>(p)++;
				phases.get(ph_q).getProp<0>(q)++;

				// Next particle
				++Np;
			}

			// Next particle
			++it;
		}
	}

	// Merge the ghost part with the real particles
	phases.get(0).ghost_put<add_,0>();
	phases.get(1).ghost_put<add_,0>();
	phases.get(2).ghost_put<add_,0>();
	phases.get(3).ghost_put<add_,0>();

	//! \cond [create sym multi-phase multi verlet] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * ## Multi-phase cell-list usage ##
	 *
	 * The multiphase cell-list that we constructed before to create a Verlet-list can be also used
	 * directly. In this case we accumulate on the property
	 * 0 of the phase 0 the distance of the near particles from all the phases
	 *
	 * \snippet Vector/4_multiphase_celllist_verlet/main.cpp cl usage
	 *
	 */

	//! \cond [cl usage] \endcond

	// we create a reference to phase0 particle for convenience
	vector_dist<3,float, aggregate<double,double> > & current_phase = phases.get(0);

	// Get the iterator of the particles of phase 0
	auto it2 = current_phase.getIterator();

	// For each particle ...
	while (it2.isNext())
	{
		// ... p
		auto p = it2.get();

		// Get the position of the particle p
		Point<3,float> xp = current_phase.getPos(p);

		// Reset to zero the propety 0
		current_phase.getProp<0>(p) = 0.0;

		// Get an iterator of all the particles neighborhood of p
		auto Np = CL_all.getNNIteratorBox(CL_all.getCell(current_phase.getPos(p)));

		// For each particle near p
		while (Np.isNext())
		{
			// Get the particle q near to p
			auto q = Np.getP();

			// Get from which phase it come from
			auto ph_q = Np.getV();

			Point<3,float> xq = phases.get(ph_q).getPos(q);

			// we accumulate all the distances
			current_phase.getProp<0>(p) = norm(xp - xq);

			++Np;
		}

		// Next particle p
		++it2;
	}

	//! \cond [cl usage] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * ## Finalize ## {#finalize}
	 *
	 *  At the very end of the program we have always de-initialize the library
	 *
	 * \snippet Vector/4_multiphase_celllist_verlet/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_4_mp_cl Vector 4 Multi Phase cell-list
	 *
	 * # Full code # {#code}
	 *
	 * \include Vector/4_multiphase_celllist_verlet/main.cpp
	 *
	 */
}




