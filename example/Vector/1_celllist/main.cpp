
#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"

/*!
 * \page Vector_1_celllist Vector 1 Cell-list and Verlet
 *
 * [TOC]
 *
 *
 * # Vector Cell-list and Verlet # {#e1_cl}
 *
 * This example show cell lists for the distributed vector and how to add particles on a distributed vector
 *
 *
 */

int main(int argc, char* argv[])
{
	/*!
	 * \page Vector_1_celllist Vector 1 Cell-list
	 *
	 * ## Initialization ##
	 *
	 * Here we Initialize the library, we create a Box that define our domain, boundary conditions, ghost
	 *
	 * \see \ref e0_s_init
	 *
	 * \snippet Vector/1_celllist/main.cpp Initialization and parameters
	 *
	 */

	//! \cond [Initialization and parameters] \endcond

	openfpm_init(&argc,&argv);
	Vcluster<> & v_cl = create_vcluster();

	// we will place the particles on a grid like way with 128 particles on each direction
	size_t sz[3] = {128,128,128};

	// The domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(1.0/(128-2));

	//! \cond [Initialization and parameters] \endcond

	/*!
	 * \page Vector_1_celllist Vector 1 Cell-list
	 *
	 * ## %Vector create ##
	 *
	 * Here we define a distributed vector in 3D, containing 3 properties, a
	 * scalar double, a vector double[3], and a tensor or rank 2 double[3][3].
	 * In this case the vector contain 0 particles initially
	 *
	 * \see \ref vector_inst
	 *
	 * \snippet Vector/1_celllist/main.cpp vector inst
	 *
	 */

	//! \cond [vector inst] \endcond

	vector_dist<3,float, aggregate<double,double[3],double[3][3]> > vd(0,box,bc,ghost);

	//! \cond [vector inst] \endcond

	/*!
	 * \page Vector_1_celllist Vector 1 Cell-list
	 *
	 * ## Grid iterator ## {#e1_cl_gr_it}
	 *
	 * We define a grid iterator, to create particles on a grid like way.
	 * An important note is that the grid iterator, iterate only on the
	 * local nodes for each processor for example suppose to have a domain like
	 * the one in figure
	 *
	 * \verbatim

	    +---------+
	    |* * *|* *|
	    |  2  |   |
	    |* * *|* *|
	    |   ---   |
	    |* *|* * *|
	    |   |     |
	    |* *|* * *|
	    |   |  1  |
	    |* *|* * *|
	    +---------+

	  \endverbatim
	 *
	 * divided in 2 processors, the processor 1 will iterate only on the points
	 * inside the portion of space marked with one. A note grid iterator follow the
	 * boundary condition specified in vector. For a perdiodic 2D 5x5 grid we have
	 *
	 * \verbatim

	    +---------+
	    * * * * * |
	    |         |
	    * * * * * |
	    |         |
	    * * * * * |
	    |         |
	    * * * * * |
	    |         |
	    *-*-*-*-*-+

	  \endverbatim
	 *
	 * Because the right border is equivalent to the left border, while for a non periodic we have the
	 * following distribution of points
	 *
	 * \verbatim

	    *-*-*-*-*
	    |       |
	    * * * * *
	    |       |
	    * * * * *
	    |       |
	    * * * * *
	    |       |
	    *-*-*-*-*

	  \endverbatim
	 *
	 * In this loop each processor will place particles on a grid
	 *
	 * \snippet Vector/1_celllist/main.cpp grid like part
	 *
	 */

	//! \cond [grid like part] \endcond

	auto it = vd.getGridIterator(sz);

	// For all the particles
	while (it.isNext())
	{
		vd.add();

		auto key = it.get();

		// The position of the particle is given by the point-index (0,0,0) ; (0,0,1) ... (127,127,127)
		// multiplied by the spacing
		vd.getLastPos()[0] = key.get(0) * it.getSpacing(0);
		vd.getLastPos()[1] = key.get(1) * it.getSpacing(1);
		vd.getLastPos()[2] = key.get(2) * it.getSpacing(2);

		// next point
		++it;
	}

	//! \cond [grid like part] \endcond

	/*!
	 * \page Vector_1_celllist Vector 1 Cell-list
	 *
	 * ## Ghost ## {#e1_part_ghost}
	 *
	 * Because we are doing a calculation that involve neighborhood particle. if a particle is near
	 * the boundary it can require particle on other processors. Ghost_get retrieve such information from
	 * the other processor. Ghost get by default synchronize the information about particle position and
	 * no properties. If we want to synchronize also properties we have to specify which one.
	 * For example with <0,1,2> here we synchronize the scalar property (0), the vector (1), and the rank 2 tensor (2)
	 *
	 * \warning Every ghost_get by default reset the status of the ghost so the information are lost
	 *
	 * \see \ref e1_gg_options
	 *
	 * \htmlonly
	 * <a href="#" onclick="hide_show('vector-video-5')" >Video 1</a>
	 * <div style="display:none" id="vector-video-5">
	 * <video id="vid5" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-2.mp4" type="video/mp4"></video>
	 * <script>video_anim('vid5',904,1592)</script>
	 * </div>
	 * <a href="#" onclick="hide_show('vector-video-4')" >Video 2</a>
	 * <div style="display:none" id="vector-video-4">
	 * <video id="vid4" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-7.mp4" type="video/mp4"></video>
	 * <script>video_anim('vid4',0,894)</script>
	 * </div>
	 * \endhtmlonly
	 *
	 * \snippet Vector/1_celllist/main.cpp ghost get
	 *
	 */

	//! \cond [ghost get] \endcond

	vd.ghost_get<0,1,2>();

	//! \cond [ghost get] \endcond

	/*!
	 * \page Vector_1_celllist Vector 1 Cell-list
	 *
	 * ## Cell-list ## {#e1_part_celllist}
	 *
	 * Here we get a Cell-list structure. A Cell-list structure is a convenient way to get the neighborhood of
	 * each particle. For each particle we get the neighborhood particles around it. Once constructed we can get
	 * an iterator over the neighborhood of the particle p
	 *
	 * In this code as demonstration we iterate over the particle.
	 *
	 * \see \ref
	 *
	 *  for each neighborhood particle of we calculate the distance. This distance is
	 * accumulated on the property 0. On the vector property (1) the distance vector is accumulated.
	 * While on the tensor the moment of inertia is calculated.
	 *
	 * \f$ s = \sum_{q = Neighborhood(p)} |p-q| \f$
	 *
	 * \f$ v_{i} = \sum_{q = Neighborhood(p)} (p_i-q_i) \f$
	 *
	 * \f$ t_{ij} = \sum_{q = Neighborhood(p)} (p_i-q_i)(p_j-q_j) \f$
	 *
	 * \snippet Vector/1_celllist/main.cpp celllist
	 *
	 */

	//! \cond [celllist] \endcond

	float r_cut = 1.0/(128-2);

	// Get cell list
	auto NN = vd.getCellList(r_cut);

	// Get the iterator
	auto it2 = vd.getDomainIterator();

	// For each particle ...
	while (it2.isNext())
	{
		// ... p
		auto p = it2.get();

		// Get the position of the particle p
		Point<3,float> xp = vd.getPos(p);

		// Get an iterator of all the particles neighborhood of p
		auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(p)));

		// For each particle near p
		while (Np.isNext())
		{
			// Get the particle q near to p
			auto q = Np.get();

			// Get the position of the particle q
			Point<3,float> xq = vd.getPos(q);

			// Calculate the distance vector between p and q
			Point<3,float> f = (xp - xq);

			// we sum the distance of all the particles
			vd.template getProp<0>(p) += f.norm();;

			// we sum the distance of all the particles
			vd.template getProp<1>(p)[0] += f.get(0);
			vd.template getProp<1>(p)[1] += f.get(1);
			vd.template getProp<1>(p)[2] += f.get(2);

			// Moment of inertia
			vd.template getProp<2>(p)[0][0] += (xp.get(0) - xq.get(0)) * (xp.get(0) - xq.get(0));
			vd.template getProp<2>(p)[0][1] += (xp.get(0) - xq.get(0)) * (xp.get(1) - xq.get(1));
			vd.template getProp<2>(p)[0][2] += (xp.get(0) - xq.get(0)) * (xp.get(2) - xq.get(2));
			vd.template getProp<2>(p)[1][0] += (xp.get(1) - xq.get(1)) * (xp.get(0) - xq.get(0));
			vd.template getProp<2>(p)[1][1] += (xp.get(1) - xq.get(1)) * (xp.get(1) - xq.get(1));
			vd.template getProp<2>(p)[1][2] += (xp.get(1) - xq.get(1)) * (xp.get(2) - xq.get(2));
			vd.template getProp<2>(p)[2][0] += (xp.get(2) - xq.get(2)) * (xp.get(0) - xq.get(0));
			vd.template getProp<2>(p)[2][1] += (xp.get(2) - xq.get(2)) * (xp.get(1) - xq.get(1));
			vd.template getProp<2>(p)[2][2] += (xp.get(2) - xq.get(2)) * (xp.get(2) - xq.get(2));

			++Np;
		}

		// Next particle p
		++it2;
	}

	//! \cond [celllist] \endcond

	/*!
	 * \page Vector_1_celllist Vector 1 Cell-list
	 *
	 * ## Verlet-list ## {#e1_part_verlet}
	 *
	 * If the particle does not move, or does not move that much we can create a verlet list
	 * for each particle, it internally use CellList to find the neighborhood but create a Verlet
	 * list it is still an expensive operation. Consider also that this datastructure can
	 * grow quickly its size is given in byte by Np * NN_average * 8
	 *
	 * * Np number of particle
	 * * NN_average number of neighborhood in average
	 *
	 * Before using it we strongly suggest to estimate the size of the data-structure
	 *
	 * \snippet Vector/1_celllist/main.cpp verletlist
	 *
	 */

	//! \cond [verletlist] \endcond

	auto verlet = vd.getVerlet(r_cut);

	auto it3 = vd.getDomainIterator();

	// For each particle i verlet.size() == Number of particles
	while (it3.isNext())
	{
		auto p = it3.get();

		// get the position of the particle i
		Point<3,float> xp = vd.getPos(p);

		// get the Neighborhood of p
		auto NNp = verlet.getNNIterator(p.getKey());

		// for each neighborhood j of particle to i
		while (NNp.isNext())
		{
			size_t q = NNp.get();

			// Get the position of the particle neighborhood if i
			Point<3,float> xq = vd.getPos(q);

			// Calculate the distance vector between p and q
			Point<3,float> f = (xp - xq);

			// we sum the distance of all the particles
			vd.template getProp<0>(p) += f.norm();;

			// we sum the distance of all the particles
			vd.template getProp<1>(p)[0] += f.get(0);
			vd.template getProp<1>(p)[1] += f.get(1);
			vd.template getProp<1>(p)[2] += f.get(2);

			// Moment of inertia
			vd.template getProp<2>(p)[0][0] += (xp.get(0) - xq.get(0)) * (xp.get(0) - xq.get(0));
			vd.template getProp<2>(p)[0][1] += (xp.get(0) - xq.get(0)) * (xp.get(1) - xq.get(1));
			vd.template getProp<2>(p)[0][2] += (xp.get(0) - xq.get(0)) * (xp.get(2) - xq.get(2));
			vd.template getProp<2>(p)[1][0] += (xp.get(1) - xq.get(1)) * (xp.get(0) - xq.get(0));
			vd.template getProp<2>(p)[1][1] += (xp.get(1) - xq.get(1)) * (xp.get(1) - xq.get(1));
			vd.template getProp<2>(p)[1][2] += (xp.get(1) - xq.get(1)) * (xp.get(2) - xq.get(2));
			vd.template getProp<2>(p)[2][0] += (xp.get(2) - xq.get(2)) * (xp.get(0) - xq.get(0));
			vd.template getProp<2>(p)[2][1] += (xp.get(2) - xq.get(2)) * (xp.get(1) - xq.get(1));
			vd.template getProp<2>(p)[2][2] += (xp.get(2) - xq.get(2)) * (xp.get(2) - xq.get(2));

			++NNp;
		}

		++it3;
	}

	//! \cond [verletlist] \endcond

	/*!
	 * \page Vector_1_celllist Vector 1 Cell-list
	 *
	 * ## Cell-list types ## {#e1_cls_types}
	 *
	 * Cell-list can have different types. The dafault that we get from getCellList is called
	 * CELL_MEMFAST. Exist other two types of Cell-list called CELL_MEMBAL and CELL_MEMMw.
	 * They all expose the same interface and same functionality so from a user prospective there
	 * is no difference. Internally they use different memory layout and they can consume
	 * more or less memory at the cost of speed. The first called CELL_MEMFAST has memory complexity
	 * equivalent to \f$N_{cells} * N_mco\f$ where \f$N_{cells}\f$ is the number of cells and \f$N_{mco}\f$
	 * is the maximum number of particles in a cell it is the fastest in access but consume a lot
	 * of memory. The CELL_MEMBAL has memory complexity \f$N_{cells} + N_{parts}\f$ it consume less memory
	 * it **can (must be evaluated case by case)** be slower because access require the resolution of
	 *  pointer to pointer. CELL_MEMMW has memory complexity \f$N_{parts}\f$ access the cell list requires
	 *  the calculation of an hashing function. The peculiarity of this Cell-list is that does not depend
	 *  from the number of cells. We suggest in anycase to use the default until memory is not a problem.
	 *
	 * \note CELL_MEMBAL and CELLMEMMW are macro that require two parameters
	 *       one is the dimensionality of the space the other is the type of the space
	 *
	 *
	 * \snippet Vector/1_celllist/main.cpp cell_list_types
	 *
	 */

	//! \cond [cell_list_types] \endcond

	// Get cell list
	auto NN4 = vd.getCellList<CELL_MEMBAL<3,float>>(r_cut);
	auto NN5 = vd.getCellList<CELL_MEMMW<3,float>>(r_cut);

	//! \cond [cell_list_types] \endcond

	/*!
	 * \page Vector_1_celllist Vector 1 Cell-list
	 *
	 * ## Finalize ## {#finalize}
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Vector/1_celllist/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_1_celllist Vector 1 Cell-list
	 *
	 * # Full code # {#code}
	 *
	 * \include Vector/1_celllist/main.cpp
	 *
	 */
}




