/*
 * ### WIKI 1 ###
 *
 * ## Simple example
 *
 * In this example we show 1D PSE derivative function approximation
 *
 * ### WIKI END ###
 *
 */

#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "PSE/Kernels.hpp"
#include "data_type/aggregate.hpp"
#include <cmath>


/*
 * ### WIKI 2 ###
 *
 * Here we define the function x*e^(-x*x) and its
 * second derivative in analytic form
 *
 * 2x*(2*x-3)*e^(-x^2)
 *
 */

double f_xex2(double x)
{
	return x*exp(-x*x);
}

double f_xex2(Point<1,double> & x)
{
	return x.get(0)*exp(-x.get(0)*x.get(0));
}

double Lapf_xex2(Point<1,double> & x)
{
	return 2.0*x.get(0)*(2.0*x.get(0)*x.get(0) - 3.0)*exp(-x.get(0)*x.get(0));
}

/*
 *
 * ### WIKI END ###
 *
 */

int main(int argc, char* argv[])
{
	//
	// ### WIKI 3 ###
	//
	// Some useful parameters. Like
	//
	// * Number of particles
	// * Minimum number of padding particles
	// * The computational domain
	// * The spacing
	// * The mollification length
	// * Second order Laplacian kernel in 1D
	//

	// Number of particles
	const size_t Npart = 125;

	// Number of padding particles (At least)
	const size_t Npad = 40;

	// The domain
	Box<1,double> box({0.0},{4.0});

	// Calculated spacing
	double spacing = box.getHigh(0) / Npart;

	// Epsilon of the particle kernel
	const double eps = 2*spacing;

	// Laplacian PSE kernel 1 dimension, on double, second order
	Lap_PSE<1,double,2> lker(eps);

	//
	// ### WIKI 2 ###
	//
	// Here we Initialize the library and we define Ghost size
	// and non-periodic boundary conditions
	//
	openfpm_init(&argc,&argv);
	Vcluster<> & v_cl = create_vcluster();

    size_t bc[1]={NON_PERIODIC};
	Ghost<1,double> g(12*eps);

	//
	// ### WIKI 3 ###
	//
	// Here we are creating a distributed vector defined by the following parameters
	//
	// we create a set of N+1 particles to have a fully covered domain of particles between 0.0 and 4.0
	// Suppose we have a spacing given by 1.0 you need 4 +1 particles to cover your domain
	//
	vector_dist<1,double, aggregate<double> > vd(Npart+1,box,bc,g);

	//
	// ### WIKI 4 ###
	//
	// We assign the position to the particles, the scalar property is set to
	// the function x*e^(-x*x) value.
	// Each processor has parts of the particles and fill part of the space, the
	// position is assigned independently from the decomposition.
	//
	// In this case, if there are 1001 particles and 3 processors the in the
	// domain from 0.0 to 4.0
	//
	// * processor 0 place particles from 0.0 to 1.332 (334 particles)
	// * processor 1 place particles from 1.336 to 2.668 (334 particles)
	// * processor 2 place particles from 2.672 to 4.0 (333 particles)
	//

	// It return how many particles the processors with id < rank has in total
	size_t base = vd.init_size_accum(Npart+1);
	auto it2 = vd.getIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		// set the position of the particles
		vd.getPos(key)[0] = (key.getKey() + base) * spacing;
		//set the property of the particles
		vd.template getProp<0>(key) = f_xex2((key.getKey() + base) * spacing);

		++it2;
	}

	//
	// ### WIKI 5 ###
	//
	// Once defined the position, we distribute them across the processors
	// following the decomposition, finally we get the ghost part
	//
	vd.map();
	vd.ghost_get<0>();

	//
	// ### WIKI 6 ###
	//
	// near the boundary we have two options, or we use one sided kernels,
	// or we add additional particles, it is required that such particles
	// produces a 2 time differentiable function. In order to obtain such
	// result we extend for x < 0.0 and x > 4.0 with the test function xe^(-x*x).
	//
	// Note that for x < 0.0 such extension is equivalent to mirror the
	// particles changing the sign of the strength
	//
	// \verbatim
	//
	// 0.6  -0.5      0.5  0.6   Strength
	//  +----+----*----*----*-
	//          0.0              Position
	//
	//  with * = real particle
	//       + = mirror particle
	//
	// \endverbatim
	//
	//
	Box<1,double> m_pad({0.0},{0.1});
	Box<1,double> m_pad2({3.9},{4.0});
	double enlarge = 0.1;

	// Create a box
	if (Npad * spacing > 0.1)
	{
		m_pad.setHigh(0,Npad * spacing);
		m_pad2.setLow(0,4.0 - Npad*spacing);
		enlarge = Npad * spacing;
	}

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		// set the position of the particles
		if (m_pad.isInsideNB(vd.getPos(key)) == true)
		{
			vd.add();
			vd.getLastPos()[0] = - vd.getPos(key)[0];
			vd.template getLastProp<0>() = - vd.template getProp<0>(key);
		}

		// set the position of the particles
		if (m_pad2.isInsideNB(vd.getPos(key)) == true)
		{
			vd.add();
			vd.getLastPos()[0] = 2.0 * box.getHigh(0) - vd.getPos(key)[0];
			vd.template getLastProp<0>() = f_xex2(vd.getLastPos()[0]);
		}

		++it;
	}

	//
	// ### WIKI 7 ###
	//
	// We create a CellList with cell spacing 12 sigma
	//

    // get and construct the Cell list
	auto cl = vd.getCellList(12*eps, CL_NON_SYMMETRIC, false, enlarge);

    // Maximum infinity norm
    double linf = 0.0;

	//
	// ### WIKI 8 ###
	//
    // For each particle get the neighborhood of each particle
    //
    // This cycle is literally the formula from PSE operator approximation
    //
    // \$ \frac{1}{\epsilon^{2}} h (u_{q} - u_{p}) \eta_{\epsilon}(x_q - x_p) \$
    //

    auto it_p = vd.getDomainIterator();
    while (it_p.isNext())
    {
    	// double PSE integration accumulator
    	double pse = 0;

    	// key
    	vect_dist_key_dx key = it_p.get();

    	// Get the position of the particles
    	Point<1,double> p = vd.getPos(key);

    	// We are not interested in calculating out the domain
    	// note added padding particle are considered domain particles
    	if (p.get(0) < 0.0 || p.get(0) >= 4.0)
    	{
    		++it_p;
    		continue;
    	}

    	// Get f(x) at the position of the particle
    	double prp_x = vd.template getProp<0>(key);

    	// Get the neighborhood of the particles
	auto NN = cl.getNNIteratorBox(cl.getCell(p));
    	while(NN.isNext())
    	{
    		auto nnp = NN.get();

    		// Calculate contribution given by the kernel value at position p,
    		// given by the Near particle
    		if (nnp != key.getKey())
    		{
    			// W(x-y)
    			double ker = lker.value(p,vd.getPos(nnp));

    			// f(y)
    			double prp_y = vd.template getProp<0>(nnp);

    			// 1.0/(eps)^2 [f(y)-f(x)]*W(x,y)*V_q
    			double prp = 1.0/eps/eps * (prp_y - prp_x) * spacing;
    			pse += prp * ker;
    		}

    		// Next particle
    		++NN;
    	}

    	// Here we calculate the L_infinity norm or the maximum difference
    	// of the analytic solution from the PSE calculated

    	double sol = Lapf_xex2(p);

    	if (fabs(pse - sol) > linf)
    		linf = fabs(pse - sol);

    	++it_p;
    }

	//
	// ### WIKI 9 ###
	//
    // Calculate the maximum infinity norm across processors and
    // print it
    //

    v_cl.max(linf);
    v_cl.execute();

    if (v_cl.getProcessUnitID() == 0)
    	std::cout << "Norm infinity: " << linf << "\n";

	//
	// ### WIKI 10 ###
	//
	// Deinitialize the library
	//
	openfpm_finalize();
}
