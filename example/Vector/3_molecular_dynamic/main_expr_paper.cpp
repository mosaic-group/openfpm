#include "Vector/vector_dist.hpp"
#include "Plot/GoogleChart.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "timer.hpp"


void Init_grid(size_t (& sz)[3],
					 vector_dist<3,double, aggregate<Point<3,double>,Point<3,double>> > & particles)
{
	auto it = particles.getGridIterator(sz);

	while (it.isNext())
	{
		particles.add();

		auto key = it.get();

		particles.getLastPos()[0] = key.get(0) * it.getSpacing(0);
		particles.getLastPos()[1] = key.get(1) * it.getSpacing(1);
		particles.getLastPos()[2] = key.get(2) * it.getSpacing(2);

		++it;
	}
}

double sigma12;
double sigma6;
double r_cut2;

constexpr int velocity_prop = 0;
constexpr int force_prop = 1;

///// Define lennard_jones interaction to be used in applyKernel_in_sim
DEFINE_INTERACTION_3D(ln_force)
	Point<3,double> r = xp - xq;
	double rn = norm2(r);

	if (rn > r_cut2)	return 0.0;

	return 24.0*(2.0 * sigma12 / (rn*rn*rn*rn*rn*rn*rn) -  sigma6 / (rn*rn*rn*rn)) * r;
END_INTERACTION

int main(int argc, char* argv[])
{
	///// Initialize the library /////
	openfpm_init(&argc,&argv);

	///// Initialize constants and initial configuration /////
	double dt = 0.0005, sigma = 0.1, r_cut = 3.0*sigma;
	sigma6 = pow(sigma,6), sigma12 = pow(sigma,12);
	double rc2 = r_cut * r_cut;

	///// Define initialization grid, simulation box, periodicity
	///// and ghost
	size_t sz[3] = {10,10,10};
	Box<3,double> box({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};
	Ghost<3,double> ghost(r_cut);

	///// Lennard_jones potential object used in applyKernel_in
	ln_force lennard_jones;

	////// Define particle and initialize on a grid
	////// the value 0 mean we start with 0 particles, Init_grid
	////// create them
	vector_dist<3,double, aggregate<Point<3,double>,Point<3,double>> > particles(0,box,bc,ghost);
	Init_grid(sz,particles);

	// Take aliases for the particle force, velocity and position
	// we use them to write math expressions
	auto force = getV<force_prop>(particles);
	auto velocity = getV<velocity_prop>(particles);
	auto position = getV<POS_PROP>(particles);

	// set the velocity to zero (a math expression)
	velocity = 0;

	///// Get the cell-list and calculate the force using lennard_jones
	///// potential
	auto NN = particles.getCellListSym(r_cut);
	force = applyKernel_in_sim(particles,NN,lennard_jones);

	// Molecula Dynamic time stepping
	for (size_t i = 0; i < 10000 ; i++)
	{
		// 1-step Verlet velocity
		// v(t + 1/2*dt) = v(t) + 1/2*force(t)*dt
		// x(t + dt) = x(t) + v(t + 1/2*dt)
		velocity = velocity + 0.5*dt*force;
		position = position + velocity*dt;

		// redistribute particles and fill ghost
		particles.map();
		particles.ghost_get<>();

		// Calculate the force at t + dt
		particles.updateCellList(NN);
		force = applyKernel_in_sim(particles,NN,lennard_jones);

		// 2-step Verlet velocity
		// v(t+dt) = v(t + 1/2*dt) + 1/2*force(t+dt)*dt
		velocity = velocity + 0.5*dt*force;
	}

	openfpm_finalize();
}

