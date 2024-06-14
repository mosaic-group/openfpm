//#define CHECKFOR_POSNAN
//#define CHECKFOR_PROPNAN

/*!
 * \page Vector_8_DEM Vector 8 Discrete element method
 *
 *
 * [TOC]
 *
 * ## Introduction {#Dem_introduction}
 *
 * In this example we show how to implement a simulation Discrete element simulation (DEM)
 * Using the Lorenz-force contact model. We will also use dynamic load balancing to keep the
 * simulation balanced during all simulation time. The simulation result is show in the
 * figures and video below
 *
 * \htmlonly
 * <a href="#" onclick="hide_show('vector-video-3')" >Simulation video 1</a><br>
 * <div style="display:none" id="vector-video-3">
 * <video id="vid3" width="1200" height="576" controls> <source src="http://openfpm.mpi-cbg.de/web/images/examples/8_DEM/DEM_velocity.mp4" type="video/mp4"></video>
 * </div>
 * \endhtmlonly
 *
 * \htmlonly
 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/8_DEM/DEM_60.png"/>
 * \endhtmlonly
 *
 * A classical model for DEM simulations of spherical granular flows is the Silbert model,
 * it includes a Herzian contact force and an elastic deformation of the grains. Each particles
 *  has radius \f$R\f$, mass \f$m\f$, polar momentum \f$I\f$ and is represented by the location of its center of mass \f$r_{i}\f$.
 * When two particles i and j collide or are in contact, the elastic contact deformation is given by:
 * \f$\delta_{ij} = 2R-r_{ij}\f$
 *
 * where \f$\vec{r}_{ij}=\vec{r}_i - \vec{r}_j\f$ is the distance vector
 * connecting particle centers and \f$r_{ij} = \|\vec{r}_{ij}\|_2\f$ its module.
 * The normal and tangential components of the relative velocity at the point of
 * contact is given by
 *
 * \f$\vec{v}_{n_{ij}} = \left(\vec{v}_{ij}\cdot\vec{n}_{ij}\right)\vec{n}_{ij}\f$
 * \f$\vec{v}_{t_{ij}} = \vec{v}_{ij}-\vec{v}_{n_{ij}}-R\left(\vec{\omega}_i + \vec{\omega}_j \right)\times \vec{n}_{ij} \f$
 *
 * whith \f$\vec{n}_{ij}=\vec{r}_{ij}/r_{ij}\f$ is the normal unit vector in direction of the distance vector,
 * \f$\vec{\omega}_i\f$ is the angular velocity of a particle and
 * \f$\vec{v}_{ij}=\vec{v}_i-\vec{v}_j\f$ the relative velocity between the two
 * particles. The evolution of the elastic tangential displacement \f$\vec{u}_{t_{ij}}\f$
 * is integrated when two particles are in contact using.
 *
 * \f$\vec{u'}_{t_{ij}} = \vec{u}_{t_{ij}} + \vec{v}_{t_{ij}} \delta t\f$
 *
 * Where $\delta t$ is the time step size. The deformation of the contacts points is
 *  stored for each particle and for each new contact point the elastic tangential displacement
 *  is initialized with \f$\vec{u}_{t_{ij}} = 0\f$. Thus for each pair of particle interacting
 *  the normal and tangential forces become:
 *
 * \f$\vec{F}_{n_{ij}}=\sqrt{\frac{\delta_{ij}}{2R}}\,\,\left(k_n\delta_{ij}\vec{n}_{ij}-\gamma_n
 * m_{\text{eff}}\vec{v}_{n_{ij}}\right) \, ,
 * \f$
 * \f$\vec{F}_{t_{ij}}=\sqrt{\frac{\delta_{ij}}{2R}}\,\,\left(-k_t\vec{u}_{t_{ij}}-\gamma_t
 * m_{\text{eff}}\vec{v}_{t_{ij}}\right) \, ,\f$
 *
 * where \f$k_{n,t}\f$ are the elastic constants in normal and tangential direction,
 * respectively, and \f$\gamma _{n,t}\f$ the corresponding viscoelastic constants. The
 * effective collision mass is given by \f$m_{\text{eff}}=\frac{m}{2}\f$. For each
 * contact point in order to enforce Coulomb's law \f$\|\vec{F}_{t_{ij}}\|_2 < \|
 * \mu\vec{F}_{n_{ij}}\|_2\f$ the tangential force is bounded
 * by the normal component force. In particular the elastic tangential displacement \f$\vec{u}_{t_{ij}}\f$ is
 * adjusted with \f$\vec{F}_{t_{ij}} \leftarrow
 * \vec{F}_{t_{ij}}\frac{\|\mu\vec{F}_{n_{ij}}\|_2}{\|\vec{F}_{t_{ij}}\|_2}
 * \f$
 *
 * This adjustment induce a truncation of the elastic displacement. The Coulomb condition
 * is equivalent to the case where two spheres slip against each other without
 * inducing additional deformations. Thus the deformation is truncated using:
 *
 * \f$\vec{u}_{t_{ij}} =
 * -\frac{1}{k_t}\left(\vec{F}_{t_{ij}} \sqrt{\frac{2R}{\delta_{ij}}} + \gamma_t
 * m_{\text{eff}}\vec{v}_{t_{ij}}\right) \, .\f$
 *
 * Considering that each particle \f$i\f$ interact with all the particles $j$ is in touch with , the total resultant force on particle $i$ is then computed by summing the contributions of all pair particles $(i,j)$. Considering that the grains are also under the effect of the gravitational field we obtain that the total force is given by
 *
 * \f$\vec{F}_i^{\text{tot}}=m \vec{g} + \sum_j
 * \left(\vec{F}_{n_{ij}}+\vec{F}_{t_{ij}}\right) \, ,
 * \f$
 *
 * where \f$\vec{g}\f$ is the acceleration due to gravity.
 * Because particles has also rotational degree of freedoms, the total torque on particle $i$ is calculated using
 *
 * \f$\vec{T}_i^{\text{tot}} = -R \sum_j \vec{n}_{ij}\times\vec{F}_{t_{ij}} \, .\f$
 *
 * \f$\vec{r}_i\f$ and angular velocities \f$\vec{\omega}_i\f$ for each particle \f$i\f$ at
 * time step \f$n+1\f$,
 * We integrate in time the equations
 * using leap frog scheme
 * with time step given by $\delta t = 10^{-6}\,$s
 *
 * \f$\vec{v}_i^{n+1} = \vec{v}_i^n + \frac{\delta t}{m}\vec{F}_i^{\text{tot}} \, ,
 * \qquad
 * \vec{r}_i^{n+1} = \vec{r}_i^n + \delta t \vec{v}_i^{n+1}
 * \f$
 * \f$\vec{\omega}_i^{n+1} = \vec{\omega}_i^n + \frac{\delta t}{I_i}\vec{T}_i^{\text{tot}} \, ,\f$
 *
 * where \f$\vec{r}_i^{n},\vec{v}_i^{n},\vec{\omega}_i^{n}\f$ denotes respectively
 *  the position, the speed and the rotational speed of the particle $i$ at time step
 *  \f$n\f$, and \f$\delta t\f$ the time step size.
 *
 * We simulate an avalanche down an inclined plane
 *
 * ## Constants {#dem_constants}
 *
 * This simulation require several constants to work, here we define mass, radius, delta time
 * elastic constants and stop time, simulation domain and ghost.
 *
 * \snippet Vector/8_DEM/main.cpp constants_def
 *
 * ## Sand particles {#sand_part}
 *
 * The definition of sand particles require:
 *
 *  * **velocity**: particle velocity
 *  * **force**: particle force
 *  * **omega**: angular velocity
 *  * **tau**: moment of the force
 *  * **cpd**: contact point deformation
 *  * **cpi**: particle index of the contact point
 *  * **ncp**: number of contact points
 *  * **tt**: sand particle or recipient particle
 *  * **gid**: global id unique for each particles
 *
 *
 * ## Initialization
 *
 * Initialize the particle in partice we create particles to form a base
 * two walls
 *
 *  \snippet Vector/8_DEM/main.cpp init sand part
 *
 * ## DEM Loop
 *
 * The DEM Loop is composition of the previous steps
 *
 * 1) integration of the force and calculated force momenta
 *
 * \snippet Vector/8_DEM/main.cpp integration pos and angvel
 *
 * 2) Update neighborhood data-structure if needed and dynamic load balance
 *
 * \snippet Vector/8_DEM/main.cpp dynamic load balancing
 *
 * 3) Calculate forces and momenta
 *
 * \snippet Vector/8_DEM/main.cpp calculate force
 *
 * The discrete element loop follow the the step described above.
 *
 */
//#define SE_CLASS3
//#define CHECKFOR_POSNAN
//#define STOP_ON_ERROR
//#define PRINT_STACKTRACE

#ifdef TEST_RUN

#else

#endif

#include "Vector/vector_dist.hpp"
#include "Draw/DrawParticles.hpp"

int main(int argc, char* argv[])
{
    // initialize the library
	openfpm_init(&argc,&argv);

	//! \cond [constants_def] \endcond

	double m = 1;
	double R = 0.06;
	double cutoff = 2 * R;
	double skin = 0.1 * cutoff;
	constexpr int max_contacts = 36;
	constexpr int max_contacts_def = max_contacts * 3;
	double dt = 0.0001;
	double Im = 2.0 * R * R * m / 5.0;
	double k_n = 7849;
	double k_t = 7849;
	double gamma_n = 3401;
	double gamma_t = 3401;
	double m_eff = 0.5;
#ifdef TEST_RUN
	double t_stop = 0.001;
#else
	double t_stop = 16.0000;
#endif
	double mu = 0.5;
	double max_vel;
	Vcluster<> & v_cl = create_vcluster();

	Box<3,double> domain({0.0,0.0,0.0},{17.04,6.0,6.6});
	Ghost<3,double> g(cutoff + skin);

	//! \cond [constants_def] \endcond

	constexpr int velocity = 0;
	constexpr int force = 1;
	constexpr int omega = 2;
	constexpr int tau = 3;
	constexpr int cpd = 4;
	constexpr int cpi = 5;
	constexpr int ncp = 6;
	constexpr int tt = 7;
	constexpr int gid = 8;

	size_t bc[3] = {NON_PERIODIC,PERIODIC,NON_PERIODIC};

	vector_dist<3,double,aggregate<Point<3,double>,Point<3,double>,Point<3,double>,Point<3,double>,
	                                  double[max_contacts_def],int[max_contacts],
									  int,int,int>> parts(0,domain,bc,g);

	//! \cond [init sand part] \endcond

	// virtual grid
	size_t sz[3] = {143,51,56};

	Box<3,double> sand_box({1.8,0.0,0.18},{8.58,5.9999,2.7});

    // we draw the initialization
	auto sand_it = DrawParticles::DrawBox(parts,sz,domain,sand_box);

	while (sand_it.isNext())
	{
		// ... add a particle ...
		parts.add();

		// ... and set it position ...
		parts.getLastPos()[0] = sand_it.get().get(0);
		parts.getLastPos()[1] = sand_it.get().get(1);
		parts.getLastPos()[2] = sand_it.get().get(2);

		parts.getLastProp<velocity>()[0] = 0.0;
		parts.getLastProp<velocity>()[1] = 0.0;
		parts.getLastProp<velocity>()[2] = 0.0;

		parts.getLastProp<omega>()[0] = 0.0;
		parts.getLastProp<omega>()[1] = 0.0;
		parts.getLastProp<omega>()[2] = 0.0;

		parts.getLastProp<tau>()[0] = 0.0;
		parts.getLastProp<tau>()[1] = 0.0;
		parts.getLastProp<tau>()[2] = 0.0;

		parts.getLastProp<force>()[0] = 0.0;
		parts.getLastProp<force>()[1] = 0.0;
		parts.getLastProp<force>()[2] = 0.0;

		parts.getLastProp<ncp>() = 0;

		parts.getLastProp<tt>() = 0;

		++sand_it;
	}

	Box<3,double> base_box({0.06,0.0,0.06},{16.98,5.9999,0.18});

    // we draw the initialization
	auto base_it = DrawParticles::DrawBox(parts,sz,domain,base_box);

	while (base_it.isNext())
	{
		// ... add a particle ...
		parts.add();

		// ... and set it position ...
		parts.getLastPos()[0] = base_it.get().get(0);
		parts.getLastPos()[1] = base_it.get().get(1);
		parts.getLastPos()[2] = base_it.get().get(2);

		parts.getLastProp<tt>() = 1;

		++base_it;
	}

	Box<3,double> wall_front({16.86,0.0,0.06},{16.98,5.9999,6.54});

    // we draw the initialization
	auto wall_f_it = DrawParticles::DrawBox(parts,sz,domain,wall_front);

	while (wall_f_it.isNext())
	{
		// ... add a particle ...
		parts.add();

		// ... and set it position ...
		parts.getLastPos()[0] = wall_f_it.get().get(0);
		parts.getLastPos()[1] = wall_f_it.get().get(1);
		parts.getLastPos()[2] = wall_f_it.get().get(2);

		parts.getLastProp<tt>() = 1;

		++wall_f_it;
	}

	Box<3,double> wall_back({0.06,0.0,0.06},{0.18,5.9999,6.54});

    // we draw the initialization
	auto wall_b_it = DrawParticles::DrawBox(parts,sz,domain,wall_back);

	while (wall_b_it.isNext())
	{
		// ... add a particle ...
		parts.add();

		// ... and set it position ...
		parts.getLastPos()[0] = wall_b_it.get().get(0);
		parts.getLastPos()[1] = wall_b_it.get().get(1);
		parts.getLastPos()[2] = wall_b_it.get().get(2);

		parts.getLastProp<tt>() = 1;

		++wall_b_it;
	}

	//! \cond [init sand part] \endcond

	parts.map();
	parts.addComputationCosts(ModelSquare());
	parts.getDecomposition().decompose();
	parts.map();

	// Fill the gid

	auto it_p = parts.getDomainIterator();
	size_t accu = parts.accum();

	while (it_p.isNext())
	{
		auto p = it_p.get();

		parts.getProp<gid>(p) = accu;

		++accu;
		++it_p;
	}

	parts.ghost_get<>();


	size_t cnt = 0;
	size_t cnt_reb = 0;
	auto nlist = parts.getVerlet<VERLET_MEMBAL<3,double,VL_NON_SYMMETRIC>>(cutoff + skin);

	double tot_sp = 0.0;

	double t = 0.0;
	while (t < t_stop)
	{
		auto pit = parts.getDomainIterator();

		max_vel = 0.0;

		//! \cond [integration pos and angvel] \endcond

		// Update
		while (pit.isNext())
		{
			auto p = pit.get();

			if (parts.getProp<tt>(p) == 1)	{++pit;continue;}

			parts.getProp<velocity>(p)[0] = parts.getProp<velocity>(p)[0] + parts.getProp<force>(p)[0]*dt;
			parts.getProp<velocity>(p)[1] = parts.getProp<velocity>(p)[1] + parts.getProp<force>(p)[1]*dt;
			parts.getProp<velocity>(p)[2] = parts.getProp<velocity>(p)[2] + parts.getProp<force>(p)[2]*dt;

			double norm2 = parts.getProp<velocity>(p)[0]*parts.getProp<velocity>(p)[0] +
						   parts.getProp<velocity>(p)[1]*parts.getProp<velocity>(p)[1] +
						   parts.getProp<velocity>(p)[2]*parts.getProp<velocity>(p)[2];
			if (max_vel < norm2)
			{max_vel = norm2;}

			parts.getPos(p)[0] = parts.getPos(p)[0] + parts.getProp<velocity>(p)[0]*dt;
			parts.getPos(p)[1] = parts.getPos(p)[1] + parts.getProp<velocity>(p)[1]*dt;
			parts.getPos(p)[2] = parts.getPos(p)[2] + parts.getProp<velocity>(p)[2]*dt;

		    if (parts.getPos(p)[0] < domain.getLow(0) || parts.getPos(p)[0] > domain.getHigh(0) ||
				parts.getPos(p)[2] < domain.getLow(2) || parts.getPos(p)[2] > domain.getHigh(2) )
		    {parts.getProp<tt>(p) = 1;}

		    parts.getProp<omega>(p)[0] = parts.getProp<omega>(p)[0] + parts.getProp<tau>(p)[0]/Im*dt;
		    parts.getProp<omega>(p)[1] = parts.getProp<omega>(p)[1] + parts.getProp<tau>(p)[1]/Im*dt;
		    parts.getProp<omega>(p)[2] = parts.getProp<omega>(p)[2] + parts.getProp<tau>(p)[2]/Im*dt;

			++pit;
		}
		tot_sp += sqrt(max_vel)*dt;
		v_cl.max(tot_sp);
		v_cl.execute();

		//! \cond [integration pos and angvel] \endcond

		//! \cond [dynamic load balancing] \endcond

		if (tot_sp >= skin / 2.0)
		{
			parts.map();

			// Check if it is time to rebalance

			if (cnt_reb >= 200)
			{
				if (v_cl.rank() == 0)
				{std::cout << "Redecomposing" << std::endl;}
				cnt_reb = 0;
				parts.addComputationCosts(ModelSquare());
				parts.getDecomposition().redecompose(200);
				parts.map();
			}

			if (v_cl.rank() == 0)
			{std::cout << "Reconstruct Verlet" << std::endl;}

			parts.ghost_get<velocity,omega,gid>();
			parts.updateVerlet(nlist,cutoff+skin);

			tot_sp = 0.0;
		}
		else
		{
			parts.ghost_get<velocity,omega,gid>();
		}

		//! \cond [dynamic load balancing] \endcond

		//! \cond [calculate force] \endcond

		auto pit2 = parts.getDomainIterator();

		while (pit2.isNext())
		{
			auto p = pit2.get();
			Point<3,double> xp = parts.getPos(p);
			Point<3,double> v_p = parts.getProp<velocity>(p);
			Point<3,double> omega_p = parts.getProp<omega>(p);

			if (parts.getProp<tt>(p) == 1)	{++pit2;continue;}

			Point<3,double> dF_n({0.0,0.0,0.0});
			Point<3,double> dF_t({0.0,0.0,0.0});
			Point<3,double> dTau({0.0,0.0,0.0});

			int contact_ok[max_contacts];
			for (size_t i = 0; i < max_contacts ; i++)	{contact_ok[i] = 0;}

			auto NN = nlist.getNNIterator(p.getKey());

			while (NN.isNext())
			{
				auto q = NN.get();

				if (q == p.getKey())	{++NN;continue;}

				Point<3,double> xq = parts.getPos(q);
				Point<3,double> v_q = parts.getProp<velocity>(q);
				Point<3,double> omega_q = parts.getProp<omega>(q);

				Point<3,double> r_pq = xp - xq;
				double r_s_pq2 = norm2(r_pq);

				// Norm is not defined, next particle
				if (r_s_pq2 == 0)
				{continue;}

				double delta_ij = 2.0*R - sqrt(r_s_pq2);

				if (delta_ij < 0.0)	{++NN;continue;}

				size_t cnt_end = parts.getProp<ncp>(p);
				int this_contact = cnt_end;

				for (size_t k = 0 ; k < cnt_end ; k++)
				{
					if (parts.getProp<cpi>(p)[k] == parts.getProp<gid>(q))
					{
						this_contact = k;
					}
				}

				int cidx;
				if (this_contact == cnt_end)
				{
					parts.getProp<ncp>(p) += 1;
					this_contact = parts.getProp<ncp>(p) - 1;

					cidx = 3 * this_contact;

					if (this_contact >= max_contacts)
					{std::cout << "Error reached maximum nunber of contacts points" << std::endl;}

					parts.getProp<cpi>(p)[this_contact] = parts.getProp<gid>(q);

					parts.getProp<cpd>(p)[cidx] = 0.0;
					parts.getProp<cpd>(p)[cidx + 1] = 0.0;
					parts.getProp<cpd>(p)[cidx + 2] = 0.0;

				}
				else
				{
					cidx = 3 * this_contact;
				}

				Point<3,double> n_ij = r_pq / sqrt(r_s_pq2);
				Point<3,double> v_rel = v_p - v_q;
				Point<3,double> v_nij = (v_rel * n_ij) * n_ij;
				Point<3,double> v_omega = (omega_p + omega_q)*R;
				Point<3,double> v_cross;

				v_cross.get(0) = v_omega.get(1) * n_ij.get(2) - v_omega.get(2) * n_ij.get(1);
				v_cross.get(1) = v_omega.get(2) * n_ij.get(0) - v_omega.get(0) * n_ij.get(2);
				v_cross.get(2) = v_omega.get(0) * n_ij.get(1) - v_omega.get(1) * n_ij.get(0);

				Point<3,double> v_tij = v_rel - v_nij - v_cross;
				Point<3,double> v_dtij = dt * v_tij;

				parts.getProp<cpd>(p)[cidx] += v_dtij.get(0);
				parts.getProp<cpd>(p)[cidx + 1] += v_dtij.get(1);
				parts.getProp<cpd>(p)[cidx + 2] += v_dtij.get(2);

				Point<3,double> u_ij;

				u_ij.get(0) = parts.getProp<cpd>(p)[cidx];
				u_ij.get(1) = parts.getProp<cpd>(p)[cidx + 1];
				u_ij.get(2) = parts.getProp<cpd>(p)[cidx + 2];

				Point<3,double> F_nij = sqrt(delta_ij/2/R) * (k_n*delta_ij*n_ij - gamma_n*m_eff*v_nij);
				dF_n = dF_n + F_nij;

				Point<3,double> F_tij = sqrt(delta_ij/2/R) * (-k_t*u_ij - gamma_t*m_eff*v_tij);
				double F_tij_sq = norm2(F_tij);
				double F_nij_sq = mu * mu * norm2(F_nij);
				if (F_tij_sq > F_nij_sq)
				{
					F_tij = F_tij * (F_nij_sq / F_tij_sq);

					parts.getProp<cpd>(p)[cidx] = -1.0 / k_t * (F_tij.get(0) * sqrt(2*R/delta_ij) + gamma_t*m_eff*v_tij.get(0));
					parts.getProp<cpd>(p)[cidx+1] = -1.0 / k_t * (F_tij.get(1) * sqrt(2*R/delta_ij) + gamma_t*m_eff*v_tij.get(1));
					parts.getProp<cpd>(p)[cidx+2] = -1.0 / k_t * (F_tij.get(2) * sqrt(2*R/delta_ij) + gamma_t*m_eff*v_tij.get(2));
				}


		        dF_t = dF_t + F_tij;
		        dTau.get(0) = dTau.get(0) - R * (n_ij.get(1) * F_tij.get(2) - n_ij.get(2) * F_tij.get(1));
		        dTau.get(1) = dTau.get(1) - R * (n_ij.get(2) * F_tij.get(0) - n_ij.get(0) * F_tij.get(2));
		        dTau.get(2) = dTau.get(2) - R * (n_ij.get(0) * F_tij.get(1) - n_ij.get(1) * F_tij.get(0));

		        contact_ok[this_contact] = 1;

				++NN;
			}

		    int cnt_end = parts.getProp<ncp>(p);
		    int i = 0;
		    for (int iread = 0; iread < cnt_end ; iread++)
		    {
		    	if (contact_ok[iread] == 1)
		    	{
		    		i = i + 1;
		    		int j = 3*(i - 1);
					int k = 3*iread;

					parts.getProp<cpd>(p)[j] = parts.getProp<cpd>(p)[k];
					parts.getProp<cpd>(p)[j+1] = parts.getProp<cpd>(p)[k+1];
					parts.getProp<cpd>(p)[j+2] = parts.getProp<cpd>(p)[k+2];
		    	}
		    }

		    parts.getProp<ncp>(p) = i;

		    if (parts.getProp<tt>(p) == 0)
		    {
		        parts.getProp<force>(p).get(0) = m * 4.905 + dF_n.get(0) + dF_t.get(0);
		        parts.getProp<force>(p).get(1) = 0.0 + dF_n.get(1) + dF_t.get(1);
		        parts.getProp<force>(p).get(2) = m * -8.49570921 + dF_n.get(2) + dF_t.get(2);

		        parts.getProp<tau>(p) = dTau;
		    }

		    if (parts.getProp<tt>(p) == 1)
		    {
		    	parts.getProp<force>(p) = 0;
		    	parts.getProp<tau>(p) = 0;
		    }

			++pit2;
		}

		//! \cond [calculate force] \endcond

		if (v_cl.rank() == 0)
		{std::cout << "Time step" << std::endl;}

		if (cnt % 300 == 0)
		{
			std::cout << "Write " << cnt << std::endl;
			parts.write_frame("output",cnt);
		}

		cnt_reb++;
		cnt++;
		t += dt;
	}

	openfpm_finalize();
}
