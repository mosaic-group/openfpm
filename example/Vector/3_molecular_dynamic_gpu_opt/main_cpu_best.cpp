#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
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


template<typename VerletList>
void calc_forces(vector_dist<3,real_number, aggregate<real_number[3],real_number[3]> > & vd, VerletList & NN, real_number sigma12, real_number sigma6, real_number r_cut)
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
	auto it2 = vd.getParticleIteratorCRS(NN);

	//! \cond [real and ghost] \endcond

	// For each particle p ...
	while (it2.isNext())
	{
		// ... get the particle p
		auto p = it2.get();

		// Get the position xp of the particle
		Point<3,real_number> xp = vd.getPos(p);

		// Get an iterator over the neighborhood particles of p
		// Note that in case of symmetric
		auto Np = NN.getNNIterator(p);

		// For each neighborhood particle ...
		while (Np.isNext())
		{
			// ... q
			auto q = Np.get();

			// if (p == q) skip this particle
			if (q == p)	{++Np; continue;};

			// Get the position of q
			Point<3,real_number> xq = vd.getPos(q);

			// Get the distance between p and q
			Point<3,real_number> r = xp - xq;

			// take the norm of this vector
			real_number rn = norm2(r);

			if (rn > r_cut * r_cut) {++Np;continue;}

			// Calculate the force, using pow is slower
			Point<3,real_number> f = 24.0*(2.0 *sigma12 / (rn*rn*rn*rn*rn*rn*rn) -  sigma6 / (rn*rn*rn*rn)) * r;

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

	// Sum the contribution to the real particles
	vd.ghost_put<add_,force>();
}


template<typename VerletList>
real_number calc_energy(vector_dist<3,real_number, aggregate<real_number[3],real_number[3]> > & vd, VerletList & NN, real_number sigma12, real_number sigma6, real_number r_cut)
{
	real_number E = 0.0;

	real_number rc = r_cut*r_cut;
	real_number shift = 4.0 * ( sigma12 / (rc*rc*rc*rc*rc*rc) - sigma6 / ( rc*rc*rc) );

	// Get an iterator over particles
	auto it2 = vd.getParticleIteratorCRS(NN);

	// For each particle ...
	while (it2.isNext())
	{
		// ... p
		auto p = it2.get();

		// Get the position of the particle p
		Point<3,real_number> xp = vd.getPos(p);

		// Get an iterator over the neighborhood particles of p
		auto Np = NN.getNNIterator(p);

		real_number Ep = E;

		// For each neighborhood of the particle p
		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			// if p == q skip this particle
			if (q == p)	{++Np; continue;};

			// Get position of the particle q
			Point<3,real_number> xq = vd.getPos(q);

			// take the normalized direction
			real_number rn = norm2(xp - xq);

			if (rn >= r_cut*r_cut)	{++Np;continue;}

			// potential energy (using pow is slower)
			E += 4.0 * ( sigma12 / (rn*rn*rn*rn*rn*rn) - sigma6 / ( rn*rn*rn) ) - shift;

			// Next neighborhood
			++Np;
		}

		// To note that the crossing scheme go across the domain particles +
		// some ghost particles. This mean that we have to filter out the ghost
		// particles otherwise we count real_number energies
		//
		if (p < vd.size_local())
		{
			// Kinetic energy of the particle given by its actual speed
			E += (vd.template getProp<velocity>(p)[0]*vd.template getProp<velocity>(p)[0] +
					vd.template getProp<velocity>(p)[1]*vd.template getProp<velocity>(p)[1] +
					vd.template getProp<velocity>(p)[2]*vd.template getProp<velocity>(p)[2]) / 2;
		}

		// Next Particle
		++it2;
	}

	// Calculated energy
	return E;
}

int main(int argc, char* argv[])
{
	real_number dt = 0.00005;
	real_number sigma = 0.01;
	real_number r_cut = 3.0*sigma;
	real_number r_gskin = 1.3*r_cut;
	real_number sigma12 = pow(sigma,12);
	real_number sigma6 = pow(sigma,6);

	openfpm::vector<real_number> x;
	openfpm::vector<openfpm::vector<real_number>> y;

	openfpm_init(&argc,&argv);
	Vcluster<> & v_cl = create_vcluster();

	// we will use it do place particles on a 10x10x10 Grid like
	size_t sz[3] = {100,100,100};

	// domain
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(r_gskin);
	ghost.setLow(0,0.0);
	ghost.setLow(1,0.0);
	ghost.setLow(2,0.0);

	vector_dist<3,real_number, aggregate<real_number[3],real_number[3]> > vd(0,box,bc,ghost,BIND_DEC_TO_GHOST);

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

	// Get the Cell list structure
	auto NN = vd.getVerletCrs(r_gskin);

	// calculate forces
	calc_forces(vd,NN,sigma12,sigma6,r_cut);
	unsigned long int f = 0;

	int cnt = 0;
	real_number max_disp = 0.0;

	// MD time stepping
	for (size_t i = 0; i < nstep ; i++)
	{
		// Get the iterator
		auto it3 = vd.getDomainIterator();

		real_number max_displ = 0.0;

		// integrate velicity and space based on the calculated forces (Step1)
		while (it3.isNext())
		{
			auto p = it3.get();

			// here we calculate v(tn + 0.5)
			vd.template getProp<velocity>(p)[0] += 0.5*dt*vd.template getProp<force>(p)[0];
			vd.template getProp<velocity>(p)[1] += 0.5*dt*vd.template getProp<force>(p)[1];
			vd.template getProp<velocity>(p)[2] += 0.5*dt*vd.template getProp<force>(p)[2];

			Point<3,real_number> disp({vd.template getProp<velocity>(p)[0]*dt,vd.template getProp<velocity>(p)[1]*dt,vd.template getProp<velocity>(p)[2]*dt});

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
			real_number energy = calc_energy(vd,NN,sigma12,sigma6,r_cut);
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

	openfpm_finalize();
}
