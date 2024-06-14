#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include "timer.hpp"

/*!
 * \page Vector_6_complex_usage Vector 6 complex usage for validation and debugging
 *
 * [TOC]
 *
 * # Complex usage for validation and debugging # {#e6_cp_deb}
 *
 * In this example we show how the flexibility of the library can be used to perform complex
 * tasks for validation and debugging. We will use a lot of feature that has been explained in
 *  the previous examples
 *
 * In a previous example we show in case of symmetric interaction between particles how to symmetric
 *  cell list could be used to speed up the calculation.
 *
 * \see not present
 *
 * In this example we validate that both computation match, in particular because the computation depend
 * from the neighborhood, we will check and validate that the neighborhood of each
 * particle is equivalent with symmetric cell list and normal cell-list.
 *
 */

int main(int argc, char* argv[])
{
	/*!
	 * \page Vector_6_complex_usage Vector 6 complex usage for validation
	 *
	 * ## Initialization ##
	 *
	 * The initialization is classical we define cutoff radius, domain, boundary conditions,
	 * ghost part and some useful constants. We also define a struct that will contain the debug information. This structure
	 * contain an id that is the global id of the particle and a point that is the the position of the
	 *  neighborhood.
	 *
	 * \see \ref e0_s_init
	 *
	 * \snippet Vector/6_complex_usage/main.cpp initialization
	 *
	 */

	//! \cond [initialization] \endcond

	openfpm_init(&argc,&argv);

	constexpr int gid = 0;
	constexpr int nn_norm = 1;
	constexpr int nn_sym = 2;

	// used to define the domain
	float L = 1000.0;

	// Domain
	Box<3,float> box({-L,-L,-L},{L,L,L});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// cut-off radius
	float r_cut = 100.0;

	// ghost
	Ghost<3,float> ghost(r_cut);

	// Point and global id
	struct point_and_gid
	{
		// global id of the particle
		size_t id;

		// Position of the neighborhood particle
		Point<3,float> xq;

		// Used to reorder the neighborhood particles by id
		bool operator<(const struct point_and_gid & pag) const
		{
			return (id < pag.id);
		}
	};

	//! \cond [initialization] \endcond

	/*!
	 * \page Vector_6_complex_usage Vector 6 complex usage for debugging
	 *
	 * ## Particle ##
	 *
	 * For this example we will use a particle with 3 properties
	 *
	 * * The first property is a global index for the particle unique across processors
	 * * The second is a vector that contain for each neighborhood the global id of the particle
	 *   and its position
	 * * The last one is again a vector equivalent to the previous but is produced with the symmetric
	 *   cell-list
	 *
	 * \snippet Vector/6_complex_usage/main.cpp particle prop
	 *
	 */

	//! \cond [particle prop] \endcond

	// Particle properties list
	typedef  aggregate<size_t,openfpm::vector<point_and_gid>,openfpm::vector<point_and_gid>> part_prop;

	//! \cond [particle prop] \endcond

	/*!
	 * \page Vector_6_complex_usage Vector 6 complex usage for debugging
	 *
	 * ## Distributed vector ##
	 *
	 * Here we create a distributed vector of 4096 particles. We initialize the particles
	 * randomly and we assign a global index to the first property of the particle
	 *
	 * The function accum return the following number
	 *
	 * \f$  \sum_{i < proc} size\_local(i) \f$
	 *
	 * where \f$ proc \f$ is the processor where we are computing the formula and \f$ size\_local(i) \f$
	 * is the number of particles the processor \f$ i \f$ has. This number is used to
	 * produce a global id of the particles
	 *
	 * \snippet Vector/6_complex_usage/main.cpp glob id part
	 *
	 */

	//! \cond [glob id part] \endcond

	// Distributed vector
	vector_dist<3,float, part_prop > vd(4096,box,bc,ghost);

	// used to calculate a global index for each particle
	size_t start = vd.accum();

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.getPos(key)[0] = 2.0*L*((float)rand()/RAND_MAX) - L;
		vd.getPos(key)[1] = 2.0*L*((float)rand()/RAND_MAX) - L;
		vd.getPos(key)[2] = 2.0*L*((float)rand()/RAND_MAX) - L;

		vd.getProp<gid>(key) = key.getKey() + start;

		++it;
	}

	//! \cond [glob id part] \endcond

	/*!
	 * \page Vector_6_complex_usage Vector 6 complex usage for debugging
	 *
	 * ## Redistribute ##
	 *
	 * Redistribute the particles and synchronize the ghosts. In this case we
	 *  are only interested in synchronizing the global id property.
	 *
	 * \snippet Vector/6_complex_usage/main.cpp map and ghost
	 *
	 */

	//! \cond [map and ghost] \endcond

	vd.map();

	// sync the ghost
	vd.ghost_get<gid>();

	//! \cond [map and ghost] \endcond

	/*!
	 * \page Vector_6_complex_usage Vector 6 complex usage for debugging
	 *
	 * ## Calculate the neighborhood of each particle ##
	 *
	 * Here we calculate the neighborhood of each particles using a simple **cell list**.
	 * In case the particle q has a distance from the particle p smaller than r_cut.
	 * We add it in the first list of neighborhood
	 *
	 * \snippet Vector/6_complex_usage/main.cpp add nn particles
	 *
	 */

	//! \cond [add nn particles] \endcond

	auto NN = vd.getCellList(r_cut);
	auto p_it = vd.getDomainIterator();

	while (p_it.isNext())
	{
		auto p = p_it.get();

		Point<3,float> xp = vd.getPos(p);

		auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(p)));

		while (Np.isNext())
		{
			auto q = Np.get();

			// skip self interaction
			if (p.getKey() == q)
			{
				++Np;
				continue;
			}

			Point<3,float> xq = vd.getPos(q);
			Point<3,float> f = (xp - xq);

			float distance = f.norm();

			// if the distance smalle than the cut-off radius add it to the neighborhood list
			if (distance < r_cut )
			{
				vd.getProp<nn_norm>(p).add();
				vd.getProp<nn_norm>(p).last().xq = xq;
				vd.getProp<nn_norm>(p).last().id = vd.getProp<0>(q);
			}

			// Next neighborhood
			++Np;
		}

		// Next particle
		++p_it;
	}

	//! \cond [add nn particles] \endcond

	/*!
	 * \page Vector_6_complex_usage Vector 6 complex usage for debugging
	 *
	 * ## Calculate the neighborhood of each particle with symmetric cell-list ##
	 *
	 * Here we calculate the neighborhood of each particle using instead a
	 * **symmetric cell list**. In case of symmetric cell-list if we find that a
	 * particle p is neighborhood of particle q we have to add p to the neighborhood
	 * of q and q to the neighborhood of p. Because q can be a ghost particle when
	 * we add the neighborhood p to q we have to transmit such information to the real
	 * owner of the particle. This can be done with the function **ghost_put**. In this
	 * case we use the operation **merge_** that add the already calculated neighborhood
	 * with the transmitted one. More in general it merge the information together.
	 *
	 * \snippet Vector/6_complex_usage/main.cpp add nn particles sym
	 *
	 */

	//! \cond [add nn particles sym] \endcond

	auto NN2 = vd.getCellListSym(r_cut);

	auto p_it2 = vd.getDomainIterator();

	while (p_it2.isNext())
	{
		auto p = p_it2.get();

		Point<3,float> xp = vd.getPos(p);

		auto Np = NN2.getNNIteratorBoxSym(NN2.getCell(vd.getPos(p)),p.getKey(),vd.getPosVector());

		while (Np.isNext())
		{
			auto q = Np.get();

			if (p.getKey() == q)
			{
				++Np;
				continue;
			}

			// repulsive

			Point<3,float> xq = vd.getPos(q);
			Point<3,float> f = (xp - xq);

			float distance = f.norm();

			// Particle should be inside r_cut range

			if (distance < r_cut )
			{
				vd.getProp<nn_sym>(p).add();
				vd.getProp<nn_sym>(q).add();

				vd.getProp<nn_sym>(p).last().xq = xq;
				vd.getProp<nn_sym>(q).last().xq = xp;
				vd.getProp<nn_sym>(p).last().id = vd.getProp<0>(q);
				vd.getProp<nn_sym>(q).last().id = vd.getProp<0>(p);
			}

			++Np;
		}

		++p_it2;
	}

	vd.ghost_put<merge_,2>();

	//! \cond [add nn particles sym] \endcond

	/*!
	 * \page Vector_6_complex_usage Vector 6 complex usage for debugging
	 *
	 * ## Cheking for validation ##
	 *
	 * Here we check that the two calculated neighborhood match. In particular,
	 * because the order of the particles does not match, we have to first reorder
	 * by global-id, and than check the list. We cannot instead easily compare the
	 * position. The reason can be seen in this figure
	 *
	 \verbatim

		  +----------------+
		  |                |
		  |           1    |
		  |           *  2 |
		  |           3  * |
		+<----- ghost *  +<-------- real
		  |* <--- real     |* <------- ghost
		  |4               |4
		  |                |
		  |                |
		  +----------------+



	 \endverbatim
	 *
	 * In the case we are calculating the neighborhood of the particle \f$ + \f$ in case of normal Cell-list.
	 * The real particle 1,2,3 are added, while the particle 4 is added with the ghost particle coordinates.
	 * In the case of symmetric cell-list the real 1,2,3 are added to \f$ + \f$. The particle 4
	 * instead is added by the ghost put. In particular 4 real add the particle to the \f$ + \f$ ghost particle
	 * and than ghost put merge the information. This mean that the symmetric add the particle 4 with the real coordinates
	 * of the particle 4
	 *
	 * \snippet Vector/6_complex_usage/main.cpp checking
	 *
	 */

	//! \cond [checking] \endcond

	auto p_it3 = vd.getDomainIterator();

	bool ret = true;
	while (p_it3.isNext())
	{
		auto p = p_it3.get();

		vd.getProp<nn_norm>(p).sort();
		vd.getProp<nn_sym>(p).sort();

		ret &= vd.getProp<nn_norm>(p).size() == vd.getProp<nn_sym>(p).size();

		for (size_t i = 0 ; i < vd.getProp<1>(p).size() ; i++)
		{
			ret &= vd.getProp<nn_norm>(p).get(i).id == vd.getProp<nn_sym>(p).get(i).id;

			if (box.isInside(vd.getProp<nn_norm>(p).get(i).xq) == true)
			{
				ret &= vd.getProp<nn_norm>(p).get(i).xq == vd.getProp<nn_sym>(p).get(i).xq;
			}
		}

		++p_it3;
	}

	if (ret != true)
	{
		std::cout << "ERROR" << std::endl;
		exit(1);
	}

	//! \cond [checking] \endcond

	/*!
	 * \page Vector_6_complex_usage Vector 6 complex usage for validation and debugging
	 *
	 * ## Finalize ##
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Vector/6_complex_usage/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_6_complex_usage Vector 6 complex usage for validation and debugging
	 *
	 * # Full code # {#code}
	 *
	 * \include Vector/6_complex_usage/main.cpp
	 *
	 */
}
