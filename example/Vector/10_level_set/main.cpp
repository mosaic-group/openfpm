#define SE_CLASS1

#include "DCPSE/DCPSE_op/DCPSE_op.hpp"
#include "Vector/vector_dist.hpp"
#include "level_set/particle_cp/particle_cp.hpp"
#include <chrono>
#include "ellipsoid_helpfunctions.h"
#include "kernel_functions.h"


// properties of distributed vector
const int sdf = 0;	// level-set function
const int normal = 1;	// surface normal
const int curvature = 2;// (mean) curvature
const int cp = 3;	// closest point on surface
const int surf_flag = 4;// surface flag that can help debugging redistancing
const int rho = 5;	// density
const int p = 6;	// pressure
const int force = 7;	// force acting on particle
const int vel = 8;	// velocity of particle
const int vel_prev = 9;	// velocity of particle at previous dt
const int pos_prev = 10;// position of particle at previous dt

typedef aggregate<double, double[2], double, double[2], int, double, double, Point<2, double>, Point<2, double>, Point<2, double>, Point<2, double>> props;
//                |        |         |         |    	|      |         |       |	       |		 |			|
//               sdf  normal curvature closest point surf_flag  density pressure force  velocity 	previous velocity	previous position
typedef vector_dist<2, double, props> particles;

const std::string filebasename = "Geometry_2Ddroplet";

// simulation parameters:

// ellipse parameters
const double A = 0.75;
const double B = 0.5;
// initial interparticle spacing
const double dp = 1/32.0;
// smoothing length of SPH operators
const double H = 3.0*dp;
// smoothing length of continuum surface
const double H_1D = H;
// polytropic exponent of equation of state p(\rho)
const double gamma_eos = 7.0;
// speed of sound
const double c = 100.0;
// reference density (of both fluids)
const double rho_0 = 1.0;
// viscosity (of both fluids)
const double eta = 0.5;
const double nu = eta*rho_0;
// surface tension coefficient
const double alpha = 50.0;
// time step size
const double dt = 2.0*1e-4;
// time variable
double t;
// variable indicating whether we are in a predictor or corrector step
int corr;

// parameters of the particle closest point redistancing method
const double band_width = 15.0*dp;
const double rdist_cutoff_factor = 2.8;

// dimensions of spatial and temporal domain
const double l = 2.44;
#ifdef TEST_RUN
const double t_end = 0.1;
#else
const double t_end = 1.0;
#endif
// total mass in the domain to compute individual particle masses from
const double M = l*l*rho_0;
// number of particles in total
double np;
// individual particle masses
double m;

// initialize kernel functions, both for the SPH operators and the continuum surface
kernel_function <2, WendlandC2>kernel2D(H);
kernel_function <1, WendlandC2>kernel1D(H_1D);

// SPH density summation to determine the individual particle densities (and volumes)
template <typename CellList> inline void density_summation(particles & vd, CellList & NN)
{
    vd.template ghost_get<rho>();
    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);
    // iterate over particles a (central particle)
    while (part.isNext())
    {
        auto a = part.get();
        Point<2, double> xa = vd.getPos(a);
        // intitialize sum that yields 1/(particle volume)
	double V_inv = 0.0;

        auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(a)));
	// iterate over particles b (neighboring particles)
        while (Np.isNext() == true)
        {
            auto b = Np.get();
            Point<2, double> xb = vd.getPos(b);
            Point<2, double> dr = xa - xb;
            double r2 = norm2(dr);
	    // if particles interact, compute and add contribution
            if (r2 <= 4.0*H*H)
            {
                double r = std::sqrt(r2);
                V_inv += kernel2D.Wab(r);
            }
            ++Np;
        }
	// particle density = (particle mass)/(particle volume)
        vd.getProp<rho>(a) = m*V_inv;
        ++part;
    }
}
// Equation of state to compute particle pressures from particle densities
// Evaluate entire list of particles in one go
inline void EqOfState(particles & vd)
{
    auto it = vd.getDomainIterator();
    // for all particles a
    while (it.isNext())
    {
        auto a = it.get();
        double dens =  vd.getProp<rho>(a);
	// Equation of state
        vd.getProp<p>(a) = (c*c*rho_0/gamma_eos)*(std::pow(dens/rho_0, gamma_eos) - 1);

        ++it;
    }
}
// force computation: pressure gradient + viscous forces + surface tension force
template<typename CellList> inline void calc_forces(particles & vd, CellList & NN)
{
    vd.template ghost_get<rho, p, vel, normal>();
    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    double max_p = 0.0;
    double max_eta = 0.0;
    double avg_alpha = 0.0;
    double maxvel = 0.0;
    int numparticles = 0;
    int numsurfparticles = 0;

    // Iterate over all central particles a
    while (part.isNext())
    {
        auto a = part.get();

	// initialize variables (forces acting on particle a, and contributions due to
	// pressure gradient and viscous forces
        vd.getProp<force>(a)[0] = 0.0;
        vd.getProp<force>(a)[1] = 0.0;

        Point<2, double> p_force;
        Point<2, double> eta_force;
        for(int k;k<2;k++) p_force[k] = 0.0;
        for(int k;k<2;k++) eta_force[k] = 0.0;

        // Get the position xp of the particle
        Point<2, double> xa = vd.getPos(a);

        // Get the density and volume of the of the particle a
        double rhoa = vd.getProp<rho>(a);
        double Va = m/rhoa;

        // Get the pressure of the particle a
        double Pa = vd.getProp<p>(a);

        // Get the Velocity of the particle a
        Point<2, double> va = vd.getProp<vel>(a);

        // Get an iterator over the neighborhood particles of p
        auto Np = NN.getNNIteratorBox(NN.getCell(vd.getPos(a)));

        // For each neighborhood particle b
        while (Np.isNext() == true)
        {
            auto b = Np.get();

            // Get the position xp of the particle
            Point<2, double> xb = vd.getPos(b);

            // if (p == q) skip this particle
            if (a.getKey() == b)	{++Np; continue;};

            Point<2, double> vb = vd.getProp<vel>(b);
            double Pb = vd.getProp<p>(b);
            double rhob = vd.getProp<rho>(b);
            double Vb = m/rhob;

            // Get the distance between p and q
            Point<2, double> dr = xa - xb;
            // take the norm of this vector
            double r2 = norm2(dr);

            // if they interact
            if (r2 < 4.0*H*H)
            {
                double r = sqrt(r2);
		// compute the difference in velocities for the viscous forces
                Point<2, double> v_rel = va - vb;
		// get derivative kernel values
                Point<2, double> DW;
                double dwdrab;
                kernel2D.DWab(dr,DW,r,false,dwdrab);
		// compute pressure gradient factor
                double factor_p = - (Va*Va + Vb*Vb)*(rhoa*Pb + rhob*Pa)/(rhoa + rhob)/m;
                // compute viscosity factor
		double factor_visc = eta*(Va*Va + Vb*Vb)*dwdrab/r/m;
		// add contributions to sum
                p_force[0] += factor_p * DW.get(0);
                p_force[1] += factor_p * DW.get(1);
                eta_force[0] += factor_visc*v_rel[0];
                eta_force[1] += factor_visc*v_rel[1];

            }

            ++Np;
        }
	// compute force on particle a as sum of pressure gradient and viscous forces
        vd.getProp<force>(a)[0] += p_force[0] + eta_force[0];
        vd.getProp<force>(a)[1] += p_force[1] + eta_force[1];

	// if particle a is close enough to the surface to experience surface tension force
	if (std::abs(vd.getProp<sdf>(a)) < (2.0*H_1D))
	{
        	double stf0;
        	double stf1;
		// evaluate 1D smoothing function that distributes surface tension effect across the interface
        	double smoothing_factor = kernel1D.Wab(std::abs(vd.getProp<sdf>(a)))/rhoa;
        	// surface tension force points inwards for particles of both phases, so we need to adjust the
		// direction of one of the surface normals
		double sign_corr = 1.0;
        	if (return_sign(vd.getProp<sdf>(a) > 0)) sign_corr = -1.0;
		// continuum surface force acting on particle = smoothing function*surface tension coefficient
		// *curvature*surface normal component/density (the latter comes from the NS-eq.)
        	stf0 = smoothing_factor*alpha*vd.getProp<curvature>(a)*vd.getProp<normal>(a)[0]*sign_corr;
        	stf1 = smoothing_factor*alpha*vd.getProp<curvature>(a)*vd.getProp<normal>(a)[1]*sign_corr;

		// add surface tension force to obtain the overall force on particle a
        	vd.getProp<force>(a)[0] += stf0;
        	vd.getProp<force>(a)[1] += stf1;

		// add curvature of particle to track the average curvature value over time (can be interesting)
        	avg_alpha += vd.getProp<curvature>(a);
        	numsurfparticles++;
    	}
        
	// get and process some info to to print during the simulation
        if (va.norm() > maxvel) maxvel = va.norm();
        if (p_force.norm() > max_p) max_p = p_force.norm();
        if (eta_force.norm() > max_eta) max_eta = eta_force.norm();
        ++numparticles;
        ++part;
    }
    avg_alpha = avg_alpha/numsurfparticles;
    if ((corr) && (((int) std::round(t/dt)%5) == 0)) std::cout<<"Time step: "<<t<<", Max p force: "<<max_p<<", Max eta force: "<<max_eta<<", Max vel: "<<maxvel<<", Average curvature: "<<avg_alpha<<", number of particles: "<<numparticles<<std::endl;
    if (numparticles == 0) throw std::invalid_argument("no particles");
}
// Predictor corrector integrator
void pred_corr_int(particles & vd, double dt, int & corr)
{
    // particle iterator
    auto part = vd.getDomainIterator();

    // iterate over all particles in the domain and integrate them
    while (part.isNext())
    {
        auto a = part.get();

        if (!corr)
        {
            // store the values at time n:
            vd.getProp<pos_prev>(a) = vd.getPos(a);
            vd.getProp<vel_prev>(a) = vd.getProp<vel>(a);
            // calculate intermediate values
            Point<2, double> dx = 0.5*dt*vd.getProp<vel>(a);
            vd.getPos(a)[0] += dx[0];
            vd.getPos(a)[1] += dx[1];
            // compute acceleration due to forces that act on particle
	    vd.getProp<vel>(a) = vd.getProp<vel>(a) + 0.5*dt*vd.getProp<force>(a);

        }
        else
        {
            // correct the accelerations and velocities
            Point<2, double> x = vd.getProp<pos_prev>(a) + dt*(vd.getProp<vel_prev>(a) + 0.5*dt*vd.getProp<force>(a));
            vd.getPos(a)[0] = x[0];
            vd.getPos(a)[1] = x[1];
            vd.getProp<vel>(a) = vd.getProp<vel_prev>(a) + dt*vd.getProp<force>(a);

        }
        ++part;
    }
    corr = !corr;

}

int main(int argc, char* argv[])
{
	np = 0;
	openfpm_init(&argc, &argv);
	// initialize the domain
	Box<2, double> domain({-l/2.0, -l/2.0}, {l/2.0, l/2.0});
	size_t sz[2] = {(size_t)(l/dp),(size_t)(l/dp)};

	// periodic boundary conditions in both x- and y-direction
	size_t bc[2] = {PERIODIC, PERIODIC};
	// ghost layers required by geometric computing
	// (finding the closest sample particle for a given particle)
	Ghost<2, double> g(3.0*band_width);

	// initialize particle list
	particles vd(0, domain, bc, g, DEC_GRAN(512));

	// give properties names for writing
	openfpm::vector<std::string> names({"sdf","normal","curvature","cp","surf_flag","rho","p","f","vel","vel_prev"});
	vd.setPropNames(names);

	auto particle_it = vd.getGridIterator(sz);

	// fill particle list using a grid iterator
	while (particle_it.isNext())
	{
		double x = -l/2.0 + dp/2.0 + particle_it.get().get(0)*dp;
		double y = -l/2.0 + dp/2.0 + particle_it.get().get(1)*dp;

		vd.add();
		vd.getLastPos()[0] = x;
		vd.getLastPos()[1] = y;
		np++;

		// initialize properties of newly added particle 
		vd.template getLastProp<surf_flag>() = 0;

		double xcp = 0.0;
		double ycp = 0.0;
		// get initial SDF value by using algorithm computing distances to ellipse
		double dist = DistancePointEllipse(A, B, abs(x), abs(y), xcp, ycp);
        	dist = -1.0*return_sign(1 - sqrt((x/A)*(x/A) + (y/B)*(y/B)))*dist;

		vd.template getLastProp<sdf>() = dist;
		vd.template getLastProp<normal>()[0] = 0.0;
		vd.template getLastProp<normal>()[1] = 0.0;
		vd.template getLastProp<curvature>() = 0.0;
		vd.template getLastProp<vel>()[0] = 0.0;
		vd.template getLastProp<vel>()[1] = 0.0;
		vd.template getLastProp<vel_prev>()[0] = 0.0;
		vd.template getLastProp<vel_prev>()[1] = 0.0;
		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<rho>() = 0.0;
		vd.template getLastProp<p>() = 0.0;

		++particle_it;
	}

    Vcluster<> & v_cl = create_vcluster();
    // Get the total number of particles
    size_t tot_part = vd.size_local();
    v_cl.sum(tot_part);
    v_cl.execute();

    // compute individual particle masses
    m = M/tot_part;
    // print the CFL condition
    std::cout<<"dt should be smaller than:\n"<<0.25*H/(c+3.0)<<"\t"<<0.125*rho_0*H*H/eta<<"\t"<<0.25*sqrt(rho_0*H*H*H/(2*M_PI*alpha))<<std::endl;

    vd.map();
    vd.ghost_get<rho>();
    // get the cell list for SPH operators (we use Wendland C2, so radius=2*smoothing length)
    auto NN = vd.getCellList(2.0*H);
    // get initial densities/volumes by computing density summation
    density_summation(vd, NN);
    // get initial pressures by evaluating the equation of state
    EqOfState(vd);

    vd.deleteGhost();
    // write initial distribution
    vd.write(filebasename + "_init_before_redistancing");

    // set redistancing parameters for the initial geometric computing
    // (to compute surface normals and curvatures)
    Redist_options rdistoptions;
    // average interparticle spacing
    rdistoptions.H = dp;
    // cutoff factor for regression neighborhood
    rdistoptions.r_cutoff_factor = rdist_cutoff_factor;
    // sampling radius of how far away sample particles can be detected
    rdistoptions.sampling_radius = 0.75*band_width;
    // tolerance for both the projections onto surface during sample particle generation
    // and Newton algorithm for solving the constrained optimization problem
    rdistoptions.tolerance = 1e-12;
    // flags on which fields to compute
    rdistoptions.compute_normals = 1;
    rdistoptions.compute_curvatures = 1;
    // flag for writing the sdf. Here, we know that the initial SDF is correct, so no need
    // for writing yet.
    rdistoptions.write_sdf = 0;
    // the closest point is not known however, so why not write it
    rdistoptions.write_cp = 1;
    // polynomial degree for minter regression
    rdistoptions.minter_poly_degree = 4;
    // lp degree of polynomials
    rdistoptions.minter_lp_degree = 1.0;
    // compute number of coefficients for internal initialization
    static constexpr unsigned int num_coeffs = minter_lp_degree_one_num_coeffs(2,4);
    // initialize the redistancing object
    particle_cp_redistancing <particles, sdf, cp, normal, curvature, num_coeffs>pcprdist_init(vd, rdistoptions);
    // start recording time to measure performance
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // run redistancing
    pcprdist_init.run_redistancing();
    // write the results (same as beforehand, but now with surface normals, curvatures and closest points)
    vd.write(filebasename + "_init");

    // change the redistancing parameters for repeated redistancing during dynamic simulation
    // now we also want to write the SDF
    rdistoptions.write_sdf = 1;
    
    // initialize parameters for the time loop
    size_t write = 0;
    size_t write_interim = 0;
    size_t rdist = 0;
    size_t largebw = 0;
    size_t it = 0;
    t = 0.0;
    corr = 0;

    // enter main time loop
    while (t <= t_end)
    {
        timer it_time;

        vd.map();
        vd.ghost_get<rho>();
        vd.updateCellList(NN);

        // Calculate particle densities and from these the pressure
	density_summation(vd, NN);
        EqOfState(vd);

        if (t != 0.0)
        {
                // with a frequency of 100, so period of 0.01, redistance with a larger sampling radius to
		// ensure detection of particles entering the narrow band
		if (largebw < t*100)
                {
                    std::cout<<"Increased sampling radius for dt."<<std::endl;
                    rdistoptions.sampling_radius = 3.0*band_width;
                }
                // initialize redistancing object with current distribution of SDF values
		particle_cp_redistancing <particles, sdf, cp, normal, curvature, num_coeffs>pcprdist(vd, rdistoptions);
                // perform redistancing
		pcprdist.run_redistancing();
		// reset the sampling radius back to the original radius if the sampling radius was increased this 
		// time step
                if (largebw < t*100)
                {
                    rdistoptions.sampling_radius = 0.75*band_width;
                    largebw++;
                }
            rdist++;
        }

        // Compute forces acting on particles
        calc_forces(vd, NN);

        // Predictor step
        it++;
        pred_corr_int(vd, dt, corr);

	if (corr == 0) t +=dt;

        if (write_interim < t*10)
        {
            vd.save(filebasename + "_interim_save");
            write_interim++;
        }

	// write with a frequency of 100
        if (write < t*100.0)//(write < t*400)
        {
            vd.deleteGhost();

            vd.write_frame(filebasename, write);
            write++;

            if (v_cl.getProcessUnitID() == 0)
            {std::cout << "TIME: " << t << std::endl;}
        }

    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time required for simulation = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    openfpm_finalize();
}
