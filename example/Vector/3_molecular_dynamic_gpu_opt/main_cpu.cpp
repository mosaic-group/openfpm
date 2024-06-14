#include "Vector/vector_dist.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include "timer.hpp"

#ifdef TEST_RUN
size_t nstep = 100;
#else
size_t nstep = 1000;
#endif

typedef float real_number;

constexpr int velocity = 0;
constexpr int force = 1;


template<typename CellList> void calc_forces(vector_dist<3,real_number, aggregate<real_number[3],real_number[3]> > & vd, CellList & NN, real_number sigma12, real_number sigma6, real_number r_cut2)
{
	vd.updateCellList(NN);

	// Get an iterator over particles
	auto it2 = vd.getDomainIterator();

	// For each particle p ...
	while (it2.isNext())
	{
		// ... get the particle p
		auto p = it2.get();

		// Get the position xp of the particle
		Point<3,real_number> xp = vd.getPos(p);

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
			Point<3,real_number> xq = vd.getPos(q);

			// Get the distance between p and q
			Point<3,real_number> r = xp - xq;

			// take the norm of this vector
			real_number rn = norm2(r);

			if (rn > r_cut2)
			{++Np; continue;};

			// Calculate the force, using pow is slower
			Point<3,real_number> f = 24.0*(2.0 *sigma12 / (rn*rn*rn*rn*rn*rn*rn) -  sigma6 / (rn*rn*rn*rn)) * r;

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

template<typename CellList> real_number calc_energy(vector_dist<3,real_number, aggregate<real_number[3],real_number[3]> > & vd, CellList & NN, real_number sigma12, real_number sigma6, real_number r_cut2)
{
	real_number rc = r_cut2;
	real_number shift = 2.0 * ( sigma12 / (rc*rc*rc*rc*rc*rc) - sigma6 / ( rc*rc*rc) );

	real_number E = 0.0;

	// Get the iterator
	auto it2 = vd.getDomainIterator();

	// For each particle ...
	while (it2.isNext())
	{
		// ... p
		auto p = it2.get();

		// Get the position of the particle p
		Point<3,real_number> xp = vd.getPos(p);

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
			Point<3,real_number> xq = vd.getPos(q);

			// take the normalized direction
			real_number rn = norm2(xp - xq);

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

	return E;
}

//! \cond [calc energy2] \endcond

int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	real_number sigma = 0.01;
	real_number r_cut = 3.0*sigma;

	// we will use it do place particles on a 10x10x10 Grid like
	size_t sz[3] = {100,100,100};

	// domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(r_cut);

	real_number dt = 0.00005;
	real_number sigma12 = pow(sigma,12);
	real_number sigma6 = pow(sigma,6);

	openfpm::vector<real_number> x;
	openfpm::vector<openfpm::vector<real_number>> y;

	vector_dist<3,real_number, aggregate<real_number[3],real_number[3]> > vd(0,box,bc,ghost);

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

	timer tsim;
	tsim.start();

	//! \cond [md steps] \endcond

	// Get the Cell list structure
	auto NN = vd.getCellList(r_cut);

	// The standard
	// auto NN = vd.getCellList(r_cut);

	// calculate forces
	calc_forces(vd,NN,sigma12,sigma6,r_cut*r_cut);
	unsigned long int f = 0;

	// MD time stepping
	for (size_t i = 0; i < nstep ; i++)
	{
		timer t_s;
		t_s.start();

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
		if (i % 500 == 0)
		{
			// We write the particle position for visualization (Without ghost)
			vd.deleteGhost();
			vd.write("particles_",f);

			// we resync the ghost
			vd.ghost_get<>();

			// We calculate the energy
			real_number energy = calc_energy(vd,NN,sigma12,sigma6,r_cut*r_cut);
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

		t_s.stop();

		if (i % 100 == 0)
		{
			// resort the particle for better cache efficency
			vd.reorder(5,reorder_opt::LINEAR);
		}
	}

	tsim.stop();
	std::cout << "Time: " << tsim.getwct() << std::endl;

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

	openfpm_finalize();
}




