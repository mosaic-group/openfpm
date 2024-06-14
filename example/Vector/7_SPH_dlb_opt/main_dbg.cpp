
//#define SE_CLASS1
//#define STOP_ON_ERROR

#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"


// A constant to indicate boundary particles
#define BOUNDARY 0

// A constant to indicate fluid particles
#define FLUID 1

size_t cnt = 0;

// initial spacing between particles dp in the formulas
const double dp = 0.0085;
// Maximum height of the fluid water
// is going to be calculated and filled later on
double h_swl = 0.0;

// c_s in the formulas (constant used to calculate the sound speed)
const double coeff_sound = 20.0;

// gamma in the formulas
const double gamma_ = 7.0;

// sqrt(3.0*dp*dp) support of the kernel
const double H = 0.0147224318643;

// Eta in the formulas
const double Eta2 = 0.01 * H*H;

// alpha in the formula
const double visco = 0.1;

// cbar in the formula (calculated later)
double cbar = 0.0;

// Mass of the fluid particles
const double MassFluid = 0.000614125;

// Mass of the boundary particles
const double MassBound = 0.000614125;

// End simulation time
const double t_end = 1.5;

// Gravity acceleration
const double gravity = 9.81;

// Reference densitu 1000Kg/m^3
const double rho_zero = 1000.0;

// Filled later require h_swl, it is b in the formulas
double B = 0.0;

// Constant used to define time integration
const double CFLnumber = 0.2;

// Minimum T
const double DtMin = 0.00001;

// Minimum Rho allowed
const double RhoMin = 700.0;

// Maximum Rho allowed
const double RhoMax = 1300.0;

// Filled in initialization
double max_fluid_height = 0.0;

// Properties

// FLUID or BOUNDARY
const size_t type = 0;

// Density
const int rho = 1;

// Density at step n-1
const int rho_prev = 2;

// Pressure
const int Pressure = 3;

// Delta rho calculated in the force calculation
const int drho = 4;

// calculated force
const int force = 5;

// velocity
const int velocity = 6;

// velocity at previous step
const int velocity_prev = 7;

struct nb_f
{
	size_t id;
	double fact;

	// Used to reorder the neighborhood particles by id
	bool operator<(const struct nb_f & pag) const
	{
		return (id < pag.id);
	}
};

// Type of the vector containing particles
typedef vector_dist<3,double,aggregate<size_t,double,  double,    double,     double,     double[3], double[3], double[3],openfpm::vector<nb_f>, openfpm::vector<nb_f>, double[3], double , size_t>> particles;
//                                       |      |        |          |            |            |         |            |
//                                       |      |        |          |            |            |         |            |
//                                     type   density   density    Pressure    delta       force     velocity    velocity
//                                                      at n-1                 density                           at n - 1

struct ModelCustom
{
	template<typename Decomposition, typename vector> inline void addComputation(Decomposition & dec, const vector & vd, size_t v, size_t p)
	{
		if (vd.template getProp<type>(p) == FLUID)
			dec.addComputationCost(v,4);
		else
			dec.addComputationCost(v,3);
	}

	template<typename Decomposition> inline void applyModel(Decomposition & dec, size_t v)
	{
		dec.setSubSubDomainComputationCost(v, dec.getSubSubDomainComputationCost(v) * dec.getSubSubDomainComputationCost(v));
	}
};

inline void EqState(particles & vd)
{
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto a = it.get();

		double rho_a = vd.template getProp<rho>(a);
		double rho_frac = rho_a / rho_zero;

		vd.template getProp<Pressure>(a) = B*( rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac - 1.0);

		++it;
	}
}



const double a2 = 1.0/M_PI/H/H/H;

inline double Wab(double r)
{
	r /= H;

	if (r < 1.0)
		return (1.0 - 3.0/2.0*r*r + 3.0/4.0*r*r*r)*a2;
	else if (r < 2.0)
		return (1.0/4.0*(2.0 - r*r)*(2.0 - r*r)*(2.0 - r*r))*a2;
	else
		return 0.0;
}



const double c1 = -3.0/M_PI/H/H/H/H;
const double d1 = 9.0/4.0/M_PI/H/H/H/H;
const double c2 = -3.0/4.0/M_PI/H/H/H/H;
const double a2_4 = 0.25*a2;
// Filled later
double W_dap = 0.0;

inline void DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print)
{
	const double qq=r/H;

	if (qq < 1.0)
	{
		double qq2 = qq * qq;
		double fac = (c1*qq + d1*qq2)/r;

		DW.get(0) = fac*dx.get(0);
		DW.get(1) = fac*dx.get(1);
		DW.get(2) = fac*dx.get(2);
	}
	else if (qq < 2.0)
	{
		double wqq = (2.0 - qq);
		double fac = c2 * wqq * wqq / r;

		DW.get(0) = fac * dx.get(0);
		DW.get(1) = fac * dx.get(1);
		DW.get(2) = fac * dx.get(2);
	}
	else
	{
		DW.get(0) = 0.0;
		DW.get(1) = 0.0;
		DW.get(2) = 0.0;
	}
}


// Tensile correction
inline double Tensile(double r, double rhoa, double rhob, double prs1, double prs2)
{
	const double qq=r/H;
	//-Cubic Spline kernel
	double wab;
	if(r>H)
	{
		double wqq1=2.0f-qq;
		double wqq2=wqq1*wqq1;

		wab=a2_4*(wqq2*wqq1);
	}
	else
	{
	    double wqq2=qq*qq;
	    double wqq3=wqq2*qq;

	    wab=a2*(1.0f-1.5f*wqq2+0.75f*wqq3);
	}

	//-Tensile correction.
	double fab=wab*W_dap;
	fab*=fab; fab*=fab; //fab=fab^4
	const double tensilp1=(prs1/(rhoa*rhoa))*(prs1>0? 0.01: -0.2);
	const double tensilp2=(prs2/(rhob*rhob))*(prs2>0? 0.01: -0.2);

	return (fab*(tensilp1+tensilp2));
}


inline double Pi(const Point<3,double> & dr, double rr2, Point<3,double> & dv, double rhoa, double rhob, double massb, double & visc)
{
	const double dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
	const double dot_rr2 = dot/(rr2+Eta2);
	visc=std::max(dot_rr2,visc);

	if(dot < 0)
	{
		const float amubar=H*dot_rr2;
		const float robar=(rhoa+rhob)*0.5f;
		const float pi_visc=(-visco*cbar*amubar/robar);

		return pi_visc;
    }
	else
		return 0.0;
}


template<typename VerletList> inline double calc_forces(particles & vd, VerletList & NN, double & max_visc, double r_gskin)
{
	auto NN2 = vd.getCellList(r_gskin);

	// Reset the ghost
    auto itg = vd.getDomainAndGhostIterator();
    while (itg.isNext())
    {
        auto p = itg.get();
        // Reset force

        if (p.getKey() < vd.size_local())
        {
			// Reset the force counter (- gravity on zeta direction)
			vd.template getProp<force>(p)[0] = 0.0;
			vd.template getProp<force>(p)[1] = 0.0;
			vd.template getProp<force>(p)[2] = -gravity;
			vd.template getProp<drho>(p) = 0.0;


			// Reset the force counter (- gravity on zeta direction)
			vd.template getProp<10>(p)[0] = 0.0;
			vd.template getProp<10>(p)[1] = 0.0;
			vd.template getProp<10>(p)[2] = -gravity;
			vd.template getProp<11>(p) = 0.0;
        }
        else
        {
			vd.getProp<force>(p)[0] = 0.0;
			vd.getProp<force>(p)[1] = 0.0;
			vd.getProp<force>(p)[2] = 0.0;

			vd.getProp<drho>(p) = 0.0;
        }


        vd.getProp<8>(p).clear();
        vd.getProp<9>(p).clear();

        ++itg;
    }

    // Get an iterator over particles
    auto part = vd.getParticleIteratorCRS(NN.getInternalCellList());

	double visc = 0;

	// For each particle ...
	while (part.isNext())
	{
		// ... a
		auto a = part.get();

		// Get the position xp of the particle
		Point<3,double> xa = vd.getPos(a);

		// Take the mass of the particle dependently if it is FLUID or BOUNDARY
		double massa = (vd.getProp<type>(a) == FLUID)?MassFluid:MassBound;

		// Get the density of the of the particle a
		double rhoa = vd.getProp<rho>(a);

		// Get the pressure of the particle a
		double Pa = vd.getProp<Pressure>(a);

		// Get the Velocity of the particle a
		Point<3,double> va = vd.getProp<velocity>(a);

		// We threat FLUID particle differently from BOUNDARY PARTICLES ...
		if (vd.getProp<type>(a) != FLUID)
		{
			// If it is a boundary particle calculate the delta rho based on equation 2
			// This require to run across the neighborhoods particles of a
			auto Np = NN.getNNIterator(a);

			// For each neighborhood particle
			while (Np.isNext() == true)
			{
				// ... q
				auto b = Np.get();

				// Get the position xp of the particle
				Point<3,double> xb = vd.getPos(b);

				// if (p == q) skip this particle
				if (a == b)	{++Np; continue;};

				// get the mass of the particle
				double massb = (vd.getProp<type>(b) == FLUID)?MassFluid:MassBound;

				// Get the velocity of the particle b
				Point<3,double> vb = vd.getProp<velocity>(b);

				// Get the pressure and density of particle b
				double Pb = vd.getProp<Pressure>(b);
				double rhob = vd.getProp<rho>(b);

				// Get the distance between p and q
				Point<3,double> dr = xa - xb;
				// take the norm of this vector
				double r2 = norm2(dr);

				// If the particles interact ...
				if (r2 < 4.0*H*H)
				{
					// ... calculate delta rho
					double r = sqrt(r2);

					Point<3,double> dv = va - vb;

					Point<3,double> DW;
					DWab(dr,DW,r,false);

					const double dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
					const double dot_rr2 = dot/(r2+Eta2);
					max_visc=std::max(dot_rr2,max_visc);

					double scal = (dv.get(0)*DW.get(0)+dv.get(1)*DW.get(1)+dv.get(2)*DW.get(2));

					vd.getProp<drho>(a) += massb*scal;
					vd.getProp<drho>(b) += massa*scal;

					struct nb_f nna;
					nna.id = b;
					nna.fact = 0.0;

					struct nb_f nnb;

					// If it is a fluid we have to calculate force
					if (vd.getProp<type>(b) == FLUID)
					{
						Point<3,double> v_rel = va - vb;
						double factor = - ((vd.getProp<Pressure>(a) + vd.getProp<Pressure>(b)) / (rhoa * rhob) + Tensile(r,rhoa,rhob,Pa,Pb) + Pi(dr,r2,v_rel,rhoa,rhob,massb,visc));

						vd.getProp<force>(b)[0] -= massa * factor * DW.get(0);
						vd.getProp<force>(b)[1] -= massa * factor * DW.get(1);
						vd.getProp<force>(b)[2] -= massa * factor * DW.get(2);

						nnb.id = vd.getProp<12>(a);
						nnb.fact = -massa * factor * DW.get(2);
					}
					else
					{
						struct nb_f nnb;
						nnb.id = vd.getProp<12>(a);
						nnb.fact = 0.0;
					}

					vd.getProp<9>(a).add(nna);
					vd.getProp<9>(b).add(nnb);

				}

				++Np;
			}

			////////////////////// NN2

			// Get an iterator over the neighborhood particles of p
			auto Np2 = NN2.getNNIteratorBox(NN2.getCell(vd.getPos(a)));

			// For each neighborhood particle
			while (Np2.isNext() == true)
			{
				// ... q
				auto b = Np2.get();

				// Get the position xp of the particle
				Point<3,double> xb = vd.getPos(b);

				// if (p == q) skip this particle
				if (a == b)	{++Np2; continue;};

				// get the mass of the particle
				double massb = (vd.getProp<type>(b) == FLUID)?MassFluid:MassBound;

				// Get the velocity of the particle b
				Point<3,double> vb = vd.getProp<velocity>(b);

				// Get the pressure and density of particle b
				double Pb = vd.getProp<Pressure>(b);
				double rhob = vd.getProp<rho>(b);

				// Get the distance between p and q
				Point<3,double> dr = xa - xb;
				// take the norm of this vector
				double r2 = norm2(dr);

				// If the particles interact ...
				if (r2 < 4.0*H*H)
				{
					// ... calculate delta rho
					double r = sqrt(r2);

					Point<3,double> dv = va - vb;

					Point<3,double> DW;
					DWab(dr,DW,r,false);

					const double dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
					const double dot_rr2 = dot/(r2+Eta2);

					vd.getProp<11>(a) += massb*(dv.get(0)*DW.get(0)+dv.get(1)*DW.get(1)+dv.get(2)*DW.get(2));


					struct nb_f nna;
					nna.id = vd.getProp<12>(b);
					nna.fact = 0.0;

					vd.getProp<8>(a).add(nna);
				}

				++Np2;
			}

		}
		else
		{
			// If it is a fluid particle calculate based on equation 1 and 2

			// Get an iterator over the neighborhood particles of p
			auto Np = NN.getNNIterator(a);

			// For each neighborhood particle
			while (Np.isNext() == true)
			{
				// ... q
				auto b = Np.get();

				// Get the position xp of the particle
				Point<3,double> xb = vd.getPos(b);

				// if (p == q) skip this particle
				if (a == b)	{++Np; continue;};

				double massb = (vd.getProp<type>(b) == FLUID)?MassFluid:MassBound;
				Point<3,double> vb = vd.getProp<velocity>(b);
				double Pb = vd.getProp<Pressure>(b);
				double rhob = vd.getProp<rho>(b);

				// Get the distance between p and q
				Point<3,double> dr = xa - xb;
				// take the norm of this vector
				double r2 = norm2(dr);

				// if they interact
				if (r2 < 4.0*H*H)
				{
					double r = sqrt(r2);

					Point<3,double> v_rel = va - vb;

					Point<3,double> DW;
					DWab(dr,DW,r,false);

					double factor = - ((vd.getProp<Pressure>(a) + vd.getProp<Pressure>(b)) / (rhoa * rhob) + Tensile(r,rhoa,rhob,Pa,Pb) + Pi(dr,r2,v_rel,rhoa,rhob,massb,visc));

					vd.getProp<force>(a)[0] += massb * factor * DW.get(0);
					vd.getProp<force>(a)[1] += massb * factor * DW.get(1);
					vd.getProp<force>(a)[2] += massb * factor * DW.get(2);

					vd.getProp<force>(b)[0] -= massa * factor * DW.get(0);
					vd.getProp<force>(b)[1] -= massa * factor * DW.get(1);
					vd.getProp<force>(b)[2] -= massa * factor * DW.get(2);

					double scal = (v_rel.get(0)*DW.get(0)+v_rel.get(1)*DW.get(1)+v_rel.get(2)*DW.get(2));

					vd.getProp<drho>(a) += massb*scal;
					vd.getProp<drho>(b) += massa*scal;


					if (vd.getProp<12>(a) == 15604 && vd.getProp<12>(b) == 15601)
					{
						std::cout << "DETECTED ERROR "  << std::endl;
					}

					if (vd.getProp<12>(b) == 15604 && vd.getProp<12>(a) == 15601)
					{
						std::cout << "DETECTED ERROR " << create_vcluster().getProcessUnitID() << "   a " << a << "  b " << b << "    " << vd.size_local_with_ghost()  << std::endl;
					}


					struct nb_f nna;
					nna.id = vd.getProp<12>(b);
					nna.fact = massb * factor * DW.get(2);

					struct nb_f nnb;
					nnb.id = vd.getProp<12>(a);
					nnb.fact = -massa * factor * DW.get(2);

					vd.getProp<9>(a).add(nna);
					vd.getProp<9>(b).add(nnb);
				}

				++Np;
			}


			///////////////////// NN2


			////////////////////// NN2

			// Get an iterator over the neighborhood particles of p
			auto Np2 = NN2.getNNIteratorBox(NN2.getCell(vd.getPos(a)));

			if (a == 0 && create_vcluster().getProcessUnitID() == 0)
			{
//				std::cout << "Particle 0" << std::endl;
			}

			// For each neighborhood particle
			while (Np2.isNext() == true)
			{
				// ... q
				auto b = Np2.get();

				// Get the position xp of the particle
				Point<3,double> xb = vd.getPos(b);

				// if (p == q) skip this particle
				if (a == b)	{++Np2; continue;};

				double massb = (vd.getProp<type>(b) == FLUID)?MassFluid:MassBound;
				Point<3,double> vb = vd.getProp<velocity>(b);
				double Pb = vd.getProp<Pressure>(b);
				double rhob = vd.getProp<rho>(b);

				// Get the distance between p and q
				Point<3,double> dr = xa - xb;
				// take the norm of this vector
				double r2 = norm2(dr);

				// if they interact
				if (r2 < 4.0*H*H)
				{
					double r = sqrt(r2);

					Point<3,double> v_rel = va - vb;

					Point<3,double> DW;
					DWab(dr,DW,r,false);

					double factor = - massb*((vd.getProp<Pressure>(a) + vd.getProp<Pressure>(b)) / (rhoa * rhob) + Tensile(r,rhoa,rhob,Pa,Pb) + Pi(dr,r2,v_rel,rhoa,rhob,massb,visc));

					vd.getProp<10>(a)[0] += factor * DW.get(0);
					vd.getProp<10>(a)[1] += factor * DW.get(1);
					vd.getProp<10>(a)[2] += factor * DW.get(2);

					vd.getProp<11>(a) += massb*(v_rel.get(0)*DW.get(0)+v_rel.get(1)*DW.get(1)+v_rel.get(2)*DW.get(2));

					struct nb_f nna;
					nna.id = vd.getProp<12>(b);
					nna.fact = factor * DW.get(2);

					vd.getProp<8>(a).add(nna);
				}

				++Np2;
			}
		}

		++part;
	}

	vd.template ghost_put<add_,drho>();
	vd.template ghost_put<add_,force>();
	vd.template ghost_put<merge_,9>();
	vd.template ghost_put<add_,10>();
	vd.template ghost_put<add_,11>();

	auto part3 = vd.getDomainIterator();

	while (part3.isNext())
	{
		auto a = part3.get();

		if (vd.getProp<force>(a)[0] != vd.getProp<10>(a)[0] ||
			vd.getProp<force>(a)[1] != vd.getProp<10>(a)[1] ||
			vd.getProp<force>(a)[2] != vd.getProp<10>(a)[2])
		{
			if (create_vcluster().getProcessUnitID() == 0)
			{
				std::cout << "LISTS " << vd.getProp<12>(a) << "    "  <<  vd.getProp<8>(a).size() << "    " << vd.getProp<9>(a).size() << std::endl;
				std::cout << "ERROR: " << vd.getProp<force>(a)[2] << "   " << vd.getProp<10>(a)[2] << "   part: " << a.getKey() << std::endl;

				// Print factors and sum


				vd.getProp<8>(a).sort();
				vd.getProp<9>(a).sort();

				double sum1 = 0.0;
				double sum2 = 0.0;

				for (size_t i = 0 ; i  < vd.getProp<8>(a).size() ; i++)
				{

					std::cout << "FACT: " << vd.getProp<8>(a).get(i).fact << "  " << vd.getProp<9>(a).get(i).fact << "   " << vd.getProp<8>(a).get(i).id << "   " << vd.getProp<9>(a).get(i).id << std::endl;

					sum1 += vd.getProp<8>(a).get(i).fact;
					sum2 += vd.getProp<9>(a).get(i).fact;
				}

				std::cout << "sum: " << sum1 << " " << sum2 << std::endl;
			}

			break;
		}

		vd.getProp<8>(a).sort();
		vd.getProp<9>(a).sort();

		if (create_vcluster().getProcessUnitID() == 0)
		{
			if (vd.getProp<8>(a).size() != vd.getProp<9>(a).size())
			{
				std::cout << "ERROR: " <<  vd.getProp<12>(a) << "    " << vd.getProp<8>(a).size() << "   "  << vd.getProp<9>(a).size()  << std::endl;

				for (size_t i = 0 ; i < vd.getProp<8>(a).size() ; i++)
				{
					std::cout << "NNNN:  " << vd.getProp<9>(a).get(i).id << "    " << vd.getProp<8>(a).get(i).id << std::endl;
				}

				break;
			}
		}

		++part3;
	}
}


void max_acceleration_and_velocity(particles & vd, double & max_acc, double & max_vel)
{
	// Calculate the maximum acceleration
	auto part = vd.getDomainIterator();

	while (part.isNext())
	{
		auto a = part.get();

		Point<3,double> acc(vd.getProp<force>(a));
		double acc2 = norm2(acc);

		Point<3,double> vel(vd.getProp<velocity>(a));
		double vel2 = norm2(vel);

		if (vel2 >= max_vel)
			max_vel = vel2;

		if (acc2 >= max_acc)
			max_acc = acc2;

		++part;
	}
	max_acc = sqrt(max_acc);
	max_vel = sqrt(max_vel);

	Vcluster & v_cl = create_vcluster();
	v_cl.max(max_acc);
	v_cl.max(max_vel);
	v_cl.execute();
}


double calc_deltaT(particles & vd, double ViscDtMax)
{
	double Maxacc = 0.0;
	double Maxvel = 0.0;
	max_acceleration_and_velocity(vd,Maxacc,Maxvel);

	//-dt1 depends on force per unit mass.
	const double dt_f = (Maxacc)?sqrt(H/Maxacc):std::numeric_limits<int>::max();

	//-dt2 combines the Courant and the viscous time-step controls.
	const double dt_cv = H/(std::max(cbar,Maxvel*10.) + H*ViscDtMax);

	//-dt new value of time step.
	double dt=double(CFLnumber)*std::min(dt_f,dt_cv);
	if(dt<double(DtMin))
		dt=double(DtMin);

	return dt;
}


openfpm::vector<size_t> to_remove;


void verlet_int(particles & vd, double dt, double & max_disp)
{
	// list of the particle to remove
	to_remove.clear();

	// particle iterator
	auto part = vd.getDomainIterator();

	double dt205 = dt*dt*0.5;
	double dt2 = dt*2.0;

	max_disp = 0;

	// For each particle ...
	while (part.isNext())
	{
		// ... a
		auto a = part.get();

		// if the particle is boundary
		if (vd.template getProp<type>(a) == BOUNDARY)
		{
			// Update rho
			double rhop = vd.template getProp<rho>(a);

			// Update only the density
	    	vd.template getProp<velocity>(a)[0] = 0.0;
	    	vd.template getProp<velocity>(a)[1] = 0.0;
	    	vd.template getProp<velocity>(a)[2] = 0.0;
	    	vd.template getProp<rho>(a) = vd.template getProp<rho_prev>(a) + dt2*vd.template getProp<drho>(a);

		    vd.template getProp<rho_prev>(a) = rhop;

			++part;
			continue;
		}

		//-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
		double dx = vd.template getProp<velocity>(a)[0]*dt + vd.template getProp<force>(a)[0]*dt205;
	    double dy = vd.template getProp<velocity>(a)[1]*dt + vd.template getProp<force>(a)[1]*dt205;
	    double dz = vd.template getProp<velocity>(a)[2]*dt + vd.template getProp<force>(a)[2]*dt205;

	    vd.getPos(a)[0] += dx;
	    vd.getPos(a)[1] += dy;
	    vd.getPos(a)[2] += dz;

	    double d2 = dx*dx + dy*dy + dz*dz;

	    max_disp = (max_disp > d2)?max_disp:d2;

	    double velX = vd.template getProp<velocity>(a)[0];
	    double velY = vd.template getProp<velocity>(a)[1];
	    double velZ = vd.template getProp<velocity>(a)[2];
	    double rhop = vd.template getProp<rho>(a);

    	vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<force>(a)[0]*dt2;
    	vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<force>(a)[1]*dt2;
    	vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<force>(a)[2]*dt2;
    	vd.template getProp<rho>(a) = vd.template getProp<rho_prev>(a) + dt2*vd.template getProp<drho>(a);

	    // Check if the particle go out of range in space and in density
	    if (vd.getPos(a)[0] <  0.000263878 || vd.getPos(a)[1] < 0.000263878 || vd.getPos(a)[2] < 0.000263878 ||
	        vd.getPos(a)[0] >  0.000263878+1.59947 || vd.getPos(a)[1] > 0.000263878+0.672972 || vd.getPos(a)[2] > 0.000263878+0.903944 ||
			vd.template getProp<rho>(a) < RhoMin || vd.template getProp<rho>(a) > RhoMax)
	    {
	                   to_remove.add(a.getKey());


	                   // if we remove something the verlet are not anymore correct and must be reconstructed
	                   max_disp = 100.0;
	    }

	    vd.template getProp<velocity_prev>(a)[0] = velX;
	    vd.template getProp<velocity_prev>(a)[1] = velY;
	    vd.template getProp<velocity_prev>(a)[2] = velZ;
	    vd.template getProp<rho_prev>(a) = rhop;

		++part;
	}

	Vcluster & v_cl = create_vcluster();
	v_cl.max(max_disp);
	v_cl.execute();

	max_disp = sqrt(max_disp);

	// remove the particles
	vd.remove(to_remove,0);

	// increment the iteration counter
	cnt++;
}

void euler_int(particles & vd, double dt, double & max_disp)
{
	// list of the particle to remove
	to_remove.clear();

	// particle iterator
	auto part = vd.getDomainIterator();

	double dt205 = dt*dt*0.5;
	double dt2 = dt*2.0;

	max_disp = 0;

	// For each particle ...
	while (part.isNext())
	{
		// ... a
		auto a = part.get();

		// if the particle is boundary
		if (vd.template getProp<type>(a) == BOUNDARY)
		{
			// Update rho
			double rhop = vd.template getProp<rho>(a);

			// Update only the density
	    	vd.template getProp<velocity>(a)[0] = 0.0;
	    	vd.template getProp<velocity>(a)[1] = 0.0;
	    	vd.template getProp<velocity>(a)[2] = 0.0;
	    	vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt*vd.template getProp<drho>(a);

		    vd.template getProp<rho_prev>(a) = rhop;

			++part;
			continue;
		}

		//-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
		double dx = vd.template getProp<velocity>(a)[0]*dt + vd.template getProp<force>(a)[0]*dt205;
	    double dy = vd.template getProp<velocity>(a)[1]*dt + vd.template getProp<force>(a)[1]*dt205;
	    double dz = vd.template getProp<velocity>(a)[2]*dt + vd.template getProp<force>(a)[2]*dt205;

	    vd.getPos(a)[0] += dx;
	    vd.getPos(a)[1] += dy;
	    vd.getPos(a)[2] += dz;

	    double d2 = dx*dx + dy*dy + dz*dz;

	    max_disp = (max_disp > d2)?max_disp:d2;

	    double velX = vd.template getProp<velocity>(a)[0];
	    double velY = vd.template getProp<velocity>(a)[1];
	    double velZ = vd.template getProp<velocity>(a)[2];
	    double rhop = vd.template getProp<rho>(a);

    	vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + vd.template getProp<force>(a)[0]*dt;
    	vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + vd.template getProp<force>(a)[1]*dt;
	   	vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + vd.template getProp<force>(a)[2]*dt;
	   	vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt*vd.template getProp<drho>(a);

	    // Check if the particle go out of range in space and in density
	    if (vd.getPos(a)[0] <  0.000263878 || vd.getPos(a)[1] < 0.000263878 || vd.getPos(a)[2] < 0.000263878 ||
	        vd.getPos(a)[0] >  0.000263878+1.59947 || vd.getPos(a)[1] > 0.000263878+0.672972 || vd.getPos(a)[2] > 0.000263878+0.903944 ||
			vd.template getProp<rho>(a) < RhoMin || vd.template getProp<rho>(a) > RhoMax)
	    {
	                   to_remove.add(a.getKey());

	                   // if we remove something the verlet are not anymore correct and must be reconstructed
	                   max_disp = 100.0;
	    }

	    vd.template getProp<velocity_prev>(a)[0] = velX;
	    vd.template getProp<velocity_prev>(a)[1] = velY;
	    vd.template getProp<velocity_prev>(a)[2] = velZ;
	    vd.template getProp<rho_prev>(a) = rhop;

		++part;
	}

	// remove the particles
	vd.remove(to_remove,0);

	Vcluster & v_cl = create_vcluster();
	v_cl.max(max_disp);
	v_cl.execute();

	max_disp = sqrt(max_disp);

	// increment the iteration counter
	cnt++;
}


int main(int argc, char* argv[])
{
    // initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<3,double> domain({-0.05,-0.05,-0.05},{1.7010,0.7065,0.5025});
	size_t sz[3] = {207,90,66};

	// Fill W_dap
	W_dap = 1.0/Wab(H/1.5);

	// Here we define the boundary conditions of our problem
    size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

    double skin = 0.25 * 2*H;
    double r_gskin = 2*H + skin;

	// extended boundary around the domain, and the processor domain
    // by the support of the cubic kernel
	Ghost<3,double> g(r_gskin);

	particles vd(0,domain,bc,g,BIND_DEC_TO_GHOST);

	// You can ignore all these dp/2.0 is a trick to reach the same initialization
	// of Dual-SPH that use a different criteria to draw particles
	Box<3,double> fluid_box({dp/2.0,dp/2.0,dp/2.0},{0.4+dp/2.0,0.67-dp/2.0,0.3+dp/2.0});

	// return an iterator to the fluid particles to add to vd
	auto fluid_it = DrawParticles::DrawBox(vd,sz,domain,fluid_box);

	// here we fill some of the constants needed by the simulation
	max_fluid_height = fluid_it.getBoxMargins().getHigh(2);
	h_swl = fluid_it.getBoxMargins().getHigh(2) - fluid_it.getBoxMargins().getLow(2);
	B = (coeff_sound)*(coeff_sound)*gravity*h_swl*rho_zero / gamma_;
	cbar = coeff_sound * sqrt(gravity * h_swl);

	// for each particle inside the fluid box ...
	while (fluid_it.isNext())
	{
		// ... add a particle ...
		vd.add();

		// ... and set it position ...
		vd.getLastPos()[0] = fluid_it.get().get(0);
		vd.getLastPos()[1] = fluid_it.get().get(1);
		vd.getLastPos()[2] = fluid_it.get().get(2);

		// and its type.
		vd.template getLastProp<type>() = FLUID;

		// We also initialize the density of the particle and the hydro-static pressure given by
		//
		// rho_zero*g*h = P
		//
		// rho_p = (P/B + 1)^(1/Gamma) * rho_zero
		//

		vd.template getLastProp<Pressure>() = rho_zero * gravity *  (max_fluid_height - fluid_it.get().get(2));

		vd.template getLastProp<rho>() = pow(vd.template getLastProp<Pressure>() / B + 1, 1.0/gamma_) * rho_zero;
		vd.template getLastProp<rho_prev>() = vd.template getLastProp<rho>();
		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<velocity_prev>()[0] = 0.0;
		vd.template getLastProp<velocity_prev>()[1] = 0.0;
		vd.template getLastProp<velocity_prev>()[2] = 0.0;

		// next fluid particle
		++fluid_it;
	}


	// Recipient
	Box<3,double> recipient1({0.0,0.0,0.0},{1.6+dp/2.0,0.67+dp/2.0,0.4+dp/2.0});
	Box<3,double> recipient2({dp,dp,dp},{1.6-dp/2.0,0.67-dp/2.0,0.4+dp/2.0});

	Box<3,double> obstacle1({0.9,0.24-dp/2.0,0.0},{1.02+dp/2.0,0.36,0.45+dp/2.0});
	Box<3,double> obstacle2({0.9+dp,0.24+dp/2.0,0.0},{1.02-dp/2.0,0.36-dp,0.45-dp/2.0});
	Box<3,double> obstacle3({0.9+dp,0.24,0.0},{1.02,0.36,0.45});

	openfpm::vector<Box<3,double>> holes;
	holes.add(recipient2);
	holes.add(obstacle1);
	auto bound_box = DrawParticles::DrawSkin(vd,sz,domain,holes,recipient1);

	while (bound_box.isNext())
	{
		vd.add();

		vd.getLastPos()[0] = bound_box.get().get(0);
		vd.getLastPos()[1] = bound_box.get().get(1);
		vd.getLastPos()[2] = bound_box.get().get(2);

		vd.template getLastProp<type>() = BOUNDARY;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<rho_prev>() = rho_zero;
		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<velocity_prev>()[0] = 0.0;
		vd.template getLastProp<velocity_prev>()[1] = 0.0;
		vd.template getLastProp<velocity_prev>()[2] = 0.0;

		++bound_box;
	}

	auto obstacle_box = DrawParticles::DrawSkin(vd,sz,domain,obstacle2,obstacle1);

	while (obstacle_box.isNext())
	{
		vd.add();

		vd.getLastPos()[0] = obstacle_box.get().get(0);
		vd.getLastPos()[1] = obstacle_box.get().get(1);
		vd.getLastPos()[2] = obstacle_box.get().get(2);

		vd.template getLastProp<type>() = BOUNDARY;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<rho_prev>() = rho_zero;
		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<velocity_prev>()[0] = 0.0;
		vd.template getLastProp<velocity_prev>()[1] = 0.0;
		vd.template getLastProp<velocity_prev>()[2] = 0.0;

		++obstacle_box;
	}

	vd.map();

	// Now that we fill the vector with particles
	ModelCustom md;

	vd.addComputationCosts(md);
	vd.getDecomposition().decompose();
	vd.map();

	size_t loc = vd.size_local();
	openfpm::vector<size_t> gt;

	Vcluster & v_cl = create_vcluster();
	v_cl.allGather(loc,gt);
	v_cl.execute();
	

	auto debug_it = vd.getDomainIterator();


	size_t tot = 0;


	for (size_t i ; i < create_vcluster().getProcessUnitID() ; i++)
		tot += gt.get(i);

	while(debug_it.isNext())
	{
		auto a = debug_it.get();

		vd.getProp<12>(a) = tot;
		tot++;

		++debug_it;
	}

	vd.ghost_get<type,rho,Pressure,velocity,12>();

	auto NN = vd.getVerletCrs(r_gskin);

	// Evolve

	size_t write = 0;
	size_t it = 0;
	size_t it_reb = 0;
	double t = 0.0;
	double tot_disp = 0.0;
	double max_disp;
	while (t <= t_end)
	{
		Vcluster & v_cl = create_vcluster();
		timer it_time;

		it_reb++;
//		if (2*tot_disp >= skin)
//		{
			vd.map();

			vd.getDecomposition().write("DecBEFORE");

			if (it_reb > 5)
			{
				ModelCustom md;
				vd.addComputationCosts(md);
				vd.getDecomposition().decompose();

				vd.map();

				it_reb = 0;

				if (v_cl.getProcessUnitID() == 0)
					std::cout << "REBALANCED " << std::endl;

			}

			vd.getDecomposition().write("DecAFTER");

			// Calculate pressure from the density
			EqState(vd);

			vd.ghost_get<type,rho,Pressure,velocity,12>();

			//vd.updateVerlet(NN,r_gskin);
			NN = vd.getVerletCrs(r_gskin);

			vd.write("Debug");

			tot_disp = 0.0;

			if (v_cl.getProcessUnitID() == 0)
				std::cout << "RECONSTRUCT Verlet " << std::endl;
//		}
//		else
//		{
//			// Calculate pressure from the density
//			EqState(vd);
//
//			vd.ghost_get<type,rho,Pressure,velocity>(SKIP_LABELLING);
//		}

		double max_visc = 0.0;

		// Calc forces
		calc_forces(vd,NN,max_visc,r_gskin);

		// Get the maximum viscosity term across processors
		v_cl.max(max_visc);
		v_cl.execute();

		// Calculate delta t integration
		double dt = calc_deltaT(vd,max_visc);

		// VerletStep or euler step
		it++;
		if (it < 40)
			verlet_int(vd,dt,max_disp);
		else
		{
			euler_int(vd,dt,max_disp);
			it = 0;
		}

		tot_disp += max_disp;
		t += dt;

		if (write < t*100)
		{
//			vd.deleteGhost();
			vd.write("Geometry",write);
//			vd.ghost_get<type,rho,Pressure,velocity>(SKIP_LABELLING);
			write++;

			if (v_cl.getProcessUnitID() == 0)
				std::cout << "TIME: " << t << "  write " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "  TOT disp: " << tot_disp << "    " << cnt << std::endl;
		}
		else
		{
			if (v_cl.getProcessUnitID() == 0)
				std::cout << "TIME: " << t << "  " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "  TOT disp: " << tot_disp << "    " << cnt << std::endl;
		}
	}

	openfpm_finalize();
}
 
