/*
 * ### WIKI 1 ###
 *
 * ## Simple example
 *
 * In this example we show The 1D PSE Diffusion equation
 *
 * in particular we integrate the following equation
 *
 * \$ \frac{u}{t} = D \laplacian u \$
 *
 * ### WIKI END ###
 *
 */

#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "PSE/Kernels.hpp"
#include "Plot/util.hpp"
#include "Plot/GoogleChart.hpp"
#include "data_type/aggregate.hpp"
#include <cmath>


/*
 * ### WIKI 2 ###
 *
 * Here we define the function x*e^(-x*x) the initial condition
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


/*
 *
 * ### WIKI 3 ###
 *
 * It calculate the Laplacian of the function in one point using the PSE
 * Approximation
 *
 * \$ \frac{1}{\epsilon^{2}} h (u_{q} - u_{p}) \eta_{\epsilon}(x_q - x_p) \$
 *
 * \param p point
 * \param key vector
 * \param vd Distributed vector
 * \param eps epsilon of the kernel
 * \param spacing between particles
 * \param cl CellList
 *
 */
template<typename CellL> double calcLap(Point<1,double> p, vect_dist_key_dx key, vector_dist<1,double, aggregate<double,double> > & vd, double eps, double spacing, CellL & cl)
{
	// Laplacian PSE kernel<1,double,2> = 1 dimension, on double, second order
	Lap_PSE<1,double,2> lker(eps);

	// double PSE integration accumulator
	double pse = 0;

	// f(x_p)
	double prp_x = vd.template getProp<0>(key);

	// Get the neighborhood of the particles
	auto NN = cl.getNNIterator(cl.getCell(p));
	while(NN.isNext())
	{
		auto nnp = NN.get();

		// Calculate contribution given by the particle q near to p
		if (nnp != key.getKey())
		{
			// W(x_p-x_q) x = position of p, y = position of q
			double ker = lker.value(p,vd.getPos(nnp));

			// f(x_q)
			double prp_y = vd.template getProp<0>(nnp);

			// 1.0/(eps)^2 [f(x_q)-f(x_p)]*W(x_p-x_q)*V_q
			double prp = 1.0/eps/eps * (prp_y - prp_x) * spacing;
			pse += prp * ker;
		}

    	// Next particle
    	++NN;
	}

	// return the calculated laplacian
	return pse;
}

/*
 *
 * ### WIKI 4 ###
 *
 * It mirror the particles at the border of the domain, in particular on the left we do the operation
 * show in figure
 *
 * \verbatim
 *
 * -0.6  -0.5      0.5  0.6   Strength
     +----+----*----*----*-
	          0.0              Position

	  with * = real particle
	       + = mirror particle

   \endverbatim
 *
 * We create particle on the negative part of the domain and we assign a strength opposite to the particle in the
 * positive part
 *
 * On the right instead we prolong the initial condition (Ideally should remain almost fixed)
 *
 * \param vd vector of particles
 * \param key particle
 * \param m_pad border box(interval) on the left
 * \param m_pad2 border box(interval) on the right
 * \param box domain
 *
 */
inline void mirror(vector_dist<1,double, aggregate<double,double> > & vd, vect_dist_key_dx & key, const Box<1,double> & m_pad, Box<1,double> & m_pad2, const Box<1,double> & box)
{
	// If the particle is inside the box (interval), it mean that is at the left border
	if (m_pad.isInsideNB(vd.getPos(key)) == true)
	{
		// Add a new particle
		vd.add();

		//Set the added particle position, symmetrically positioned on the negative part of the domain
		vd.getLastPos()[0] = - vd.getPos(key)[0];

		// Set the property of the particle to be negative
		vd.template getLastProp<0>() = - vd.template getProp<0>(key);
	}

	// If the particle is inside the box (interval), it mean that is at the right border
	if (m_pad2.isInsideNB(vd.getPos(key)) == true)
	{
		// Add a new particle
		vd.add();

		//Set the added particle position, symmetrically positioned around x = 4.0
		//         4.0    Added
		//   *------|------*
		//
		vd.getLastPos()[0] = 2.0 * box.getHigh(0) - vd.getPos(key)[0];

		// no flow boundary condition
		vd.template getLastProp<0>() = vd.template getProp<0>(key);
	}
}

/*
 *
 * ### WIKI END ###
 *
 */

int main(int argc, char* argv[])
{
	//
	// ### WIKI 5 ###
	//
	// Some useful parameters.
	//

	// Number of particles
	const size_t Npart = 1000;
        
        // A different way to write Npart
        size_t NpartA[1] = {Npart+1};

	// Number of steps
	const size_t Nstep = 1000;

	// Delta t
	const double dt = 10.0 / Nstep;

	// Number of padding particles, particles outside the domain
	const size_t Npad = 20;

	// The domain
	Box<1,double> box({0.0},{4.0});

	// Calculated spacing
	double spacing = box.getHigh(0) / Npart;

	// The mollification length
	const double eps = 2*spacing;

        // Diffusion constant
	double D = 1e-4;

        // 2 constants
        constexpr int U = 0;
        constexpr int Lap_U = 1;
        
	//
	// ### WIKI 6 ###
	//
	// Here we Initialize the library and we define Ghost size
	// and non-periodic boundary conditions
	//
	openfpm_init(&argc,&argv);
	Vcluster<> & v_cl = create_vcluster();

        size_t bc[1]={NON_PERIODIC};
	Ghost<1,double> g(12*eps);

	//
	// ### WIKI 7 ###
	//
	// Here we are creating a distributed vector defined by the following parameters 
        //
        // <1,double, aggregate<double,double> > = 1 dimension, space with type double, 2 properties double,double (aggregate is like tuples, a way to encapsulate types). The first property is the field value u, the second is the Laplacian of u
	// with Npart+1 particles, on the box domain [0.0,4.0], with bc boundary condition NON_PERIODIC, and an extended ghost part of size g
	//
	vector_dist<1,double, aggregate<double,double> > vd(Npart+1,box,bc,g);

	//
	// ### WIKI 8 ###
	//
	// We assign the position to the particles in order to be equally spaced in a grid like position, the scalar property is set to
	// the function x*e^(-x*x) value the initial condition
        //
        // P.S. In a Multi processor context getGridIterator adapt based on the underline space division across processors
	//
        
        // Create a grid like iterator
        auto itp = vd.getGridIterator(NpartA);

        // Counter starting from 0
        size_t p = 0;
        
        // For each point in the grid
	while (itp.isNext())
	{
                // Get the grid index
		auto key = itp.get();
                
                // Assign the particle position
		vd.getPos(p)[0] = key.get(0) * spacing;
                
                // Assign the property 0 based on the initial condition
                vd.template getProp<U>(p) = f_xex2(key.get(0) * spacing);

                // Next particle
                ++p;
                
                // Next grid point
        	++itp;
	}

	//
	// ### WIKI 9 ###
	//
	// Once defined the position, we distribute them across the processors
	// following the decomposition, finally we update the ghost part for the vield U
	//
	vd.map();
	vd.ghost_get<U>();

	//
	// ### WIKI 10 ###
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
	//-0.6  -0.5      0.5  0.6   Strength
	//  +----+----*----*----*-
	//          0.0              Position
	//
	//  with * = real particle
	//       + = mirror particle
	//
	// \endverbatim
	//
	//
        
        // The particle that we have to mirror are ...
        
        // ... on the box m_pad for the laft part
	Box<1,double> m_pad({0.0},{0.1});
        
        // ... on the box m_pad2 for the right part
	Box<1,double> m_pad2({3.9},{4.0});
        
        // It contain how much we enlarge the domain
	double enlarge = m_pad.getHigh(0) - m_pad.getLow(0);

	// In general 0.1 should be a sufficent padding, but in case it is not we recalculate the left and right boxes
	if (Npad * spacing > 0.1)
	{
		m_pad.setHigh(0,Npad * spacing);
		m_pad2.setLow(0,4.0 - Npad*spacing);
		enlarge = Npad * spacing;
	}


	//
	// ### WIKI 11 ###
	//
	// Here we mirror the particles if needed
	//
	//
	
	auto it = vd.getDomainIterator();
	size_t n_part = vd.size_local();

	while (it.isNext())
	{
		auto key = it.get();

		mirror(vd,key,m_pad,m_pad2,box);

		++it;
	}

	//
	// ### WIKI 12 ###
	//
	// We create a CellList with cell spacing 8 epsilon, on the part of space domain + ghost + padding
	// The padding area, can be bigger than the ghost part. In practice there is no reason to do it, but
	// keeping ghost size and padding area unrelated give us the possibility to show how to create a CellList
	// on a area bigger than the domain + ghost
	
	// Create a Cell list with Cell spaping 8*epsilon, the CellList is created on a domain+ghost+enlarge space
	auto cl = vd.getCellList(8*eps, CL_NON_SYMMETRIC, false, enlarge)

	// Maximum infinity norm
	double linf = 0.0;

	//
	// ### WIKI 13 ###
	//
        // Euler time integration until t = 10.0
        //
        //

        for (double t = 0; t <= 10.0 ; t += dt)
        {
        
            // Iterate all the particles, including padding because are real particles
            auto it_p = vd.getDomainIterator();
            while (it_p.isNext())
            {
                // Get particle p
    		auto p = it_p.get();

                // Get the position of the particle p
    		Point<1,double> x_p = vd.getPos(p);

    		// We are not interested in calculating out the domain (So on the padded particles)
    		if (x_p.get(0) < 0.0 || x_p.get(0) >= 4.0)
    		{
    			++it_p;
    			continue;
    		}

    		// Here we calculate the Laplacian of u on position x_p and we store the result on the particle p
    		vd.template getProp<Lap_U>(p) = calcLap(x_p,p,vd,eps,spacing,cl);

                // Next particle
    		++it_p;
            }

            // Once we calculated the Laplacian we do not need enymore the padding particle eliminate them
            // n_part is the original size of the vector without padding particles
            vd.resize(n_part);

            // Iterate again on all particles
            auto it_p2 = vd.getDomainIterator();
            while (it_p2.isNext())
            {
                // Get particle p
    		auto p = it_p2.get();

                // Get the Laplacian
    		double pse = vd.template getProp<Lap_U>(p);

    		// Integrate with Euler step
    		vd.template getProp<0>(p) += D * pse * dt;

                // mirror again the particles if needed (this for the next step)
    		mirror(vd,p,m_pad,m_pad2,box);

                // Next particle
    		++it_p2;
            }

            // Update the ghost part
            vd.ghost_get<U>();
        }

        //
        // ### WIKI 14 ###
        //
        // U now contain the solution scattered across processors. Because it is not possible plot on multiple
        // processors (Google chart does not support it) we collect the solution on one processor
        //

        struct pos_val
        {
            // position of the particle
            double pos;
        
            // value of U of the particle
            double val;

            // usefull to sort the particle by position
            bool operator<(const pos_val & p) const
            {
                    return pos < p.pos;
            }
        };

        // Resize a vector with the local size of the vector
        openfpm::vector<pos_val> v(vd.size_local());
    
        // This will contail all the particles on the master processor
        openfpm::vector<pos_val> v_tot;
    

        // Copy the particle position and the field U, into the vector v
        for (size_t i = 0 ; i < vd.size_local() ; i++)
        {
            // particle position
            v.get(i).pos = vd.getPos(i)[0];
        
            // Field U
            v.get(i).val = vd.template getProp<U>(i);
        }

        // from the partial solution vector v we Gather the total solution v_tot on master (processor=0)
        v_cl.SGather(v,v_tot,0);

        //
        // ### WIKI 15 ###
        //
        // Once we gather the solution on master we plot it. Only processor 0 do this
        //

        // Only processor = 0
        if (v_cl.getProcessUnitID() == 0)
        {
            // sort the solution by particle position
            v_tot.sort();
        
            // y contain 2 vectors the first is the calculated solution
            // the second contain the initial condition
            openfpm::vector<openfpm::vector<double>> y(2);

            // vector with the particle position for plotting
            openfpm::vector<double> x;

            // We fill the plotting variables
            for (size_t i = 0 ; i < v_tot.size() ; i++)
            {
                // fill x with the particle position
                x.add(v_tot.get(i).pos);
                
                // fill y with the value of the field U contained by the particles
    		y.get(0).add(v_tot.get(i).val);
                
                // fill y with the initial condition
    		y.get(1).add(f_xex2(v_tot.get(i).pos));
            }

            // Plot with GooglePlot
            // Google charts options
            GCoptions options;

            // Assign a name for the Y axis
            options.yAxis = std::string("U Axis");
        
            // Assign a name for the X Axis
            options.xAxis = std::string("X Axis");
        
            // width of the line
            options.lineWidth = 1.0;

            // Create a google plot object
            GoogleChart cg;
        
            // add a line plot
            cg.AddLinesGraphT(x,y,options);
        
            // write the plot file
            cg.write("PSE_plot.html");
        }

        //
        // ### WIKI 16 ###
        //
        // Deinitialize the library
        //
        openfpm_finalize();
}
