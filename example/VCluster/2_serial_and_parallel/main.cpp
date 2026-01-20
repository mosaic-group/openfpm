#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "VCluster/VCluster.hpp"

/*!
 *
 * \page VCluster_2_serial_and_parallel Parallel and serial basics
 *
 * [TOC]
 *
 * # Parallel and serial {#vcl_2_ser_par}
 *
 * In this example we discuss how to switch from Parallel to serial and from serial to
 * parallel. This is in general useful if you are trying to port a serial program into
 * openfpm, and you want to do it progressively or in parts
 *
 */

#define POS_ 0

//! \cond [my_particle_def] \endcond

struct my_particle
{
	//! position
	double x[2];

	//! property 0
	double prop0;

	//! property 1
	long int prop1;
};

//! \cond [my_particle_def] \endcond

//! \cond [serial_to_parallel] \endcond

void serial_to_parallel(std::vector<my_particle> & parts,
						vector_dist<2,double,aggregate<double,long int>> & pparts)
{
	auto & v_cl = create_vcluster();
	pparts.clear();

	if (v_cl.getProcessUnitID() == 0)
	{
		for (size_t i = 0 ; i < parts.size() ; i++)
		{
			pparts.add();

			pparts.getLastPos()[0] = parts[i].x[0];
			pparts.getLastPos()[1] = parts[i].x[1];

			pparts.getLastProp<0>() = parts[i].prop0;
			pparts.getLastProp<1>() = parts[i].prop1;
		}
	}

	pparts.map<Error>();
	pparts.addComputationCosts();
	pparts.getDecomposition().decompose();
	pparts.map<Error>();
}

//! \cond [serial_to_parallel] \endcond

//! \cond [parallel_to_serial] \endcond

void parallel_to_serial(std::vector<my_particle> & parts,
						vector_dist<2,double,aggregate<double,long int>> & pparts)
{
	auto & v_cl = create_vcluster();

	// here we collect on the processor 0

	auto & pos = pparts.getPosVector();
	auto & prp = pparts.getPropVector();

	std::remove_reference<decltype(pos)>::type pos_collect;
	std::remove_reference<decltype(prp)>::type prp_collect;

	// we collecto everything on processor 0 on pos_collect and prp_collect
	v_cl.SGather(pos,pos_collect,0);
	v_cl.SGather(prp,prp_collect,0);

	if (v_cl.getProcessUnitID() == 0)
	{
		parts.clear();

		for (size_t i = 0 ; i < pos_collect.size() ; i++)
		{
			struct my_particle p;

			p.x[0] = pos_collect.get<POS_>(i)[0];
			p.x[1] = pos_collect.get<POS_>(i)[1];

			p.prop0 = prp_collect.get<0>(i);
			p.prop1 = prp_collect.get<1>(i);

			parts.push_back(p);
		}
	}
}

//! \cond [parallel_to_serial] \endcond

int main(int argc, char* argv[])
{
	/*!
	 *
	 * \page VCluster_2_serial_and_parallel Parallel and serial basics
	 *
	 *
	 * ## Initialization {#vcl_2_init}
	 *
	 * Before using any functionality the library must be initialized. After initialization we can create
	 * the Vcluster object
	 *
	 * \snippet VCluster/2_serial_and_parallel/main.cpp initialization
	 *
	 */

	//! \cond [initialization] \endcond

	openfpm_init(&argc,&argv);
	Vcluster<> & v_cl = create_vcluster();
	
	//! \cond [initialization] \endcond

	/*!
	 *
	 * \page VCluster_2_serial_and_parallel Parallel and serial basics
	 *
	 * ## Serial section {#vcl_2_ser_sect}
	 *
	 * First we start from a serial section in which we have one external **std::vector<my_particle>**
	 * containing particles position and properties. Because the code is serial the full serial part
	 * is enclosed into an if condition that only processor 0 execute the code. Processor 0 only execute
	 * the code and create 100 particles.
	 *
	 * \snippet VCluster/2_serial_and_parallel/main.cpp my_particle_def
	 *
	 * \snippet VCluster/2_serial_and_parallel/main.cpp serial
	 *
	 */

	//! \cond [serial] \endcond

	std::vector<my_particle> parts;

	// id of the processor calling this function
	long int proc_id = v_cl.getProcessUnitID();

	if (proc_id == 0)
	{
		// We create 100 particles randomly between (0,0) and (2.0,3.0) serially

		for (size_t i = 0 ; i < 100 ; i++)
		{
			my_particle p;

			p.x[0] = 2.0*(double)rand()/RAND_MAX;
			p.x[1] = 3.0*(double)rand()/RAND_MAX;
			p.prop0 = 0;
			p.prop1 = 0;
			parts.push_back(p);
		}
	}

	//! \cond [serial] \endcond

	/*!
	 *
	 * \page VCluster_2_serial_and_parallel Parallel and serial basics
	 *
	 * ## Parallel section {#vcl_2_par_sect}
	 *
	 * The parallel section instead require an equivalent vector_dist. In particular we
	 * initialize the vector dist with 0 particles on a domain from (0.0,0.0) to
	 * (2.0,3.0).
	 *
	 * \snippet VCluster/2_serial_and_parallel/main.cpp vector_definition
	 *
	 * As soon as we want to go parallel we want to convert the information
	 * from the **std::vector to vector_dist**. To do this, we call the function **serial_to_parallel**
	 *  that convert the serial std::vector into vector_dist. After we did the conversion
	 *  we can use **pparts** to do operations on it.
	 *
	 * ### serial_to_parallel function
	 *
	 * \snippet VCluster/2_serial_and_parallel/main.cpp serial_to_parallel
	 *
	 * This function convert the serial set of particles into the parallel set of particles
	 * To do this we first clear **pparts** than only the processor 0 fill **pparts** from the serial
	 * set of particles *parts*. After we filled **pparts** we redistribute the particles across
	 * processors using the function **map**. If the particles are uniformly distributed we can skip the
	 * second part an comment **addComputationCosts()**,**decompose()**,and the second **map()**.
	 * Otherwise these three lines redistribute the particles across processors in order to have (almost) an
	 * equal number of particles across processors independently from their density
	 *
	 * ## Serial section (again) {#vcl_2_sect_to_par_}
	 *
	 * To reconvert to serial we use the function **parallel_to_serial** this function
	 * does the opposite to convert the vector_dist **pparts** into std::vector **parts**.
	 * Once we have converted we can use **parts**, in this case we simply move the particles
	 * of 0.01 on x and 0.01 on y
	 *
	 * ### parallel_to_serial function
	 *
	 * \snippet VCluster/2_serial_and_parallel/main.cpp parallel_to_serial
	 *
	 * This function convert the parallel set of particles into the serial set of particles.
	 * To do this we first get the position and properties particle vector from each processor.
	 *  Than we create two other equivalent vector that will contain the collected vector
	 *   of position and vector of properties. We use the function **SGather** to collect
	 * on processor 0 all the data. Once collected we move the data from the collected vectors
	 *  into **parts**
	 *
	 */

	//! \cond [vector_definition] \endcond

	Box<2,double> domain({0.0,0.0},{2.0,3.0});

	Ghost<2,double> gs(0.02);
	size_t bc[2] = {PERIODIC,PERIODIC};

	vector_dist<2,double,aggregate<double,long int>> pparts(0,domain,bc,gs);

	//! \cond [vector_definition] \endcond

	// We do 100 iteration
	for (size_t i = 0 ; i < 100 ; i++)
	{
		//! \cond [serial_to_parallel] \endcond

		// here we convert from serial to parallel moving the data from
		// the serial part to pparts

		serial_to_parallel(parts,pparts);

		// Now we can use pparts to

		auto it = pparts.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			pparts.getPos(key)[0] += 0.01;
			pparts.getPos(key)[1] += 0.01;

			pparts.getProp<0>(key) = i*i;
			pparts.getProp<1>(key) = i;

			++it;
		}

		pparts.write_frame("output",i);

		//! \cond [serial_to_parallel] \endcond

		//! \cond [parallel_to_serial] \endcond

		parallel_to_serial(parts,pparts);

		// id of the processor calling this function
		long int proc_id = v_cl.getProcessUnitID();

		if (proc_id == 0)
		{
			// Here we are serial with parts updates

			for (size_t i = 0 ; i < parts.size() ; i++)
			{
				parts[i].x[0] += 0.01;
				parts[i].x[1] += 0.01;

				parts[i].prop0 = i*i;
				parts[i].prop1 = i;
			}
		}

		//! \cond [parallel_to_serial] \endcond
	}

	/*!
	 * \page VCluster_2_serial_and_parallel Parallel and serial basics
	 *
	 * ## Finalize ##
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet VCluster/1_semantic/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond
}
