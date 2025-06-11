#include "Vector/vector_dist.hpp"
#include "Plot/GoogleChart.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"

constexpr int velocity = 0;
constexpr int force = 1;

struct ln_potential
{
	double sigma12,sigma6;

	ln_potential(double sigma12_, double sigma6_) {sigma12 = sigma12_, sigma6 = sigma6_;}

	Point<2,double> value(const Point<3,double> & xp, const Point<3,double> xq)
	{
		double rn = norm2(xp - xq);

		Point<2,double> E({4.0 * ( sigma12 / (rn*rn*rn*rn*rn*rn) - sigma6 / ( rn*rn*rn) ),
						 0.0});

		return E;
	}
};

struct ln_force
{
	double sigma12,sigma6;

	ln_force(double sigma12_, double sigma6_) {sigma12 = sigma12_; sigma6 = sigma6_;}

	Point<3,double> value(const Point<3,double> & xp, const Point<3,double> xq)
	{
		Point<3,double> r = xp - xq;
		double rn = norm2(r);

		return 24.0*(2.0 * sigma12 / (rn*rn*rn*rn*rn*rn*rn) -  sigma6 / (rn*rn*rn*rn)) * r;
	}
};

int main(int argc, char* argv[])
{
	double dt = 0.0005, sigma = 0.1, r_cut = 0.3;

	double sigma6 = pow(sigma,6), sigma12 = pow(sigma,12);

	openfpm::vector<double> x;
	openfpm::vector<openfpm::vector<double>> y;


	openfpm_init(&argc,&argv);
	Vcluster & vcl = create_vcluster();

	size_t sz[3] = {10,10,10};
	Box<3,double> box({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};
	Ghost<3,double> ghost(r_cut);

	ln_force lf(sigma12,sigma6);
	ln_potential lp(sigma12,sigma6);

	vector_dist<3,double, aggregate<Point<3,double>,Point<3,double>> > vd(0,box,bc,ghost);

	auto v_force = getV<force>(vd);
	auto v_velocity = getV<velocity>(vd);
	auto v_pos = getV<POS_PROP>(vd);

	auto it = vd.getGridIterator(sz);

	while (it.isNext())
	{
		vd.add();

		auto key = it.get();

		vd.getLastPos()[0] = key.get(0) * it.getSpacing(0);
		vd.getLastPos()[1] = key.get(1) * it.getSpacing(1);
		vd.getLastPos()[2] = key.get(2) * it.getSpacing(2);

		++it;
	}

	v_force = 0;
	v_velocity = 0;

	auto NN = vd.getCellList(r_cut);

	vd.updateCellList(NN);
	v_force = applyKernel_in_sim(vd,NN,lf);
	unsigned long int f = 0;

	// MD time stepping
	for (size_t i = 0; i < 10000 ; i++)
	{
		assign(v_velocity, v_velocity + 0.5*dt*v_force,
			   v_pos, v_pos + v_velocity*dt);

		vd.map();
		vd.template ghost_get<>();

		// calculate forces or a(tn + 1) Step 2
		vd.updateCellList(NN);
		v_force = applyKernel_in_sim(vd,NN,lf);

		v_velocity = v_velocity + 0.5*dt*v_force;

		if (i % 100 == 0)
		{
			vd.deleteGhost();
			vd.write("particles_",f);
			vd.ghost_get<>();

			Point<2,double> E = rsum(applyKernel_in_sim(vd,NN,lp) + (v_velocity * v_velocity)/2.0,vd).get();

			vcl.sum(E.get(0));vcl.sum(E.get(1));
			vcl.execute();

			// we save the energy calculated at time step i c contain the time-step y contain the energy
			x.add(i);
			y.add({E.get(0),E.get(1),E.get(0) - E.get(1)});

			if (vcl.getProcessUnitID() == 0)
				std::cout << "Energy Total: " << E.get(0) << "   Kinetic: " <<  E.get(1) << "   Potential: " << E.get(0) - E.get(1) << std::endl;

			f++;
		}
	}

	GCoptions options;
	options.title = std::string("Energy with time");
	options.yAxis = std::string("Energy");
	options.xAxis = std::string("iteration");
	options.lineWidth = 1.0;

	GoogleChart cg;
	cg.AddLinesGraph(x,y,options);
	cg.write("gc_plot2_out.html");

	openfpm_finalize();
}
