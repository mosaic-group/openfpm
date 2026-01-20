/**
 * @file Vortex_in_cell/main_vic_petsc.cpp
 * @page Vortex_in_cell Vortex in cell
 *
 * @subpage Vortex_in_cell_petsc
 * @subpage Vortex_in_cell_petsc_opt
 *
 *
 */

/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Vortex in Cell 3D ring with ringlets # {#vic_ringlets}
 *
 * In this example we solve the Navier-Stokes equation in the vortex formulation in 3D
 * for an incompressible fluid. (bold symbols are vectorial quantity)
 *
 * \htmlonly
 * <a href="#" onclick="hide_show('vector-video-1')" >Video 1</a>
 * <div style="display:none" id="vector-video-1">
 * <video id="vid1" width="1200" height="576" controls> <source src="http://openfpm.mpi-cbg.de/web/images/examples/3_Vortex_in_cell/vortex_in_cell.mp4"></video>
 * <script>video_anim('vid1',100,230)</script>
 * </div>
 * <a href="#" onclick="hide_show('vector-video-2')" >Video 2</a>
 * <div style="display:none" id="vector-video-2">
 * <video id="vid2" width="1200" height="576" controls> <source src="http://openfpm.mpi-cbg.de/web/images/examples/3_Vortex_in_cell/vortex_in_cell_iso.mp4"></video>
 * <script>video_anim('vid2',21,1590)</script>
 * </div>
 * \endhtmlonly
 *
 * ## Numerical method ## {#num_vic_mt}
 *
 * In this code we solve the Navier-stokes equation for incompressible fluid in the
 * vorticity formulation. We first recall the Navier-stokes equation in vorticity formulation
 *
 * \f$ \nabla \times \boldsymbol u = - \boldsymbol w \f$
 *
 * \f$  \frac{\displaystyle D \boldsymbol w}{\displaystyle dt} = ( \boldsymbol w \cdot \vec \nabla) \boldsymbol u + \nu \nabla^{2} \boldsymbol w \f$    (5)
 *
 * Where \f$w\f$ is the vorticity and \f$u\f$ is the velocity of the fluid.
 * With Reynold number defined as \f$Re = \frac{uL}{\nu}\f$. The algorithm can be expressed with the following pseudo code.
 *
 * \verbatim

	1) Initialize the vortex ring on grid
	2) Do an helmotz hodge projection to make the vorticity divergent free
	3) Initialize particles on the same position as the grid or remesh

	while (t < t_end) do
		4) Interpolate vorticity from the particles to mesh
		5) calculate velocity u from the vorticity w
		6) calculate the right-hand-side on grid and interpolate on particles
		7) interpolate velocity u to particles
		8) move particles accordingly to the velocity
		9) interpolate the vorticity into mesh and reinitialize the particles
		   in a grid like position
	end while

   \endverbatim
 *
 * This pseudo code show how to solve the equation above using euler integration.
 * In case of Runge-kutta of order two the pseudo code change into
 *
 *
 * \verbatim

	1) Initialize the vortex ring on grid
	2) Do an helmotz hodge projection to make the vorticity divergent free
	3) Initialize particles on the same position as the grid or remesh

	while (t < t_end) do
		4) Interpolate vorticity from the particles to mesh
		5) calculate velocity u from the vorticity w
		6) calculate the right-hand-side on grid and interpolate on particles
		7) interpolate velocity u to particles
		8) move particles accordingly to the velocity and save the old position in x_old

		9) Interpolate vorticity on mesh on the particles
		10) calculate velocity u from the vorticity w
		11) calculate the right-hand-side on grid and interpolate on particles
		12) interpolate velocity u to particles
		13) move particles accordingly to the velocity starting from x_old
		14) interpolate the vorticity into mesh and reinitialize the particles
		   in a grid like position
	end while

   \endverbatim
 *
 * In the following we explain how each step is implemented in the code
 *
 * ## Inclusion ## {#num_vic_inc}
 *
 * This example code need several components. First because is a particle
 * mesh example we have to activate **grid_dist_id.hpp** and **vector_dist_id.hpp**.
 * Because we use a finite-difference scheme and linear-algebra to calculate the
 * velocity out of the vorticity, we have to include **FDScheme.hpp** to produce
 * from the finite difference scheme a matrix that represent the linear-system
 * to solve. **SparseMatrix.hpp** is the Sparse-Matrix that will contain the linear
 * system to solve in order to get the velocity out of the vorticity.
 * **Vector.hpp** is the data-structure that contain the solution of the
 * linear system. **petsc_solver.hpp** is the library to use in order invert the linear system.
 * Because we have to interpolate between particles and grid we the to include
 * **interpolate.hpp** as interpolation kernel we use the mp4, so we include the
 * **mp4_kernel.hpp**
 *
 * For convenience we also define the particles type and the grid type and some
 * convenient constants
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp inclusion
 *
 *
 *
 */

#include "config.h"

#ifdef HAVE_PETSC

//#define SE_CLASS1
//#define PRINT_STACKTRACE

//! \cond [inclusion] \endcond

#include "interpolation/interpolation.hpp"
#include "Grid/grid_dist_id.hpp"
#include "Vector/vector_dist.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Vector/Vector.hpp"
#include "FiniteDifference/FDScheme.hpp"
#include "Solvers/petsc_solver.hpp"
#include "interpolation/mp4_kernel.hpp"
#include "Solvers/petsc_solver_AMG_report.hpp"

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

// The type of the grids
typedef grid_dist_id<3,float,aggregate<float[3]>> grid_type;

// The type of the particles
typedef vector_dist<3,float,aggregate<float[3],float[3],float[3],float[3],float[3]>> particles_type;

// radius of the torus
float ringr1 = 1.0;
// radius of the core of the torus
float sigma = 1.0/3.523;
// Reynold number (If you want to use 7500.0 you have to use a grid 1600x400x400)
//float tgtre  = 7500.0;
float tgtre = 3000.0;

// Kinematic viscosity
float nu = 1.0/tgtre;

// Time step
// float dt = 0.0025 for Re 7500
float dt = 0.0125;

#ifdef TEST_RUN
const unsigned int nsteps = 10;
#else
const unsigned int nsteps = 10001;
#endif

// All the properties by index
constexpr unsigned int vorticity = 0;
constexpr unsigned int velocity = 0;
constexpr unsigned int p_vel = 1;
constexpr int rhs_part = 2;
constexpr unsigned int old_vort = 3;
constexpr unsigned int old_pos = 4;

//! \cond [inclusion] \endcond

template<typename grid> void calc_and_print_max_div_and_int(grid & g_vort)
{
	//! \cond [sanity_int_div] \endcond

	g_vort.template ghost_get<vorticity>();
	auto it5 = g_vort.getDomainIterator();

	double max_vort = 0.0;

	double int_vort[3] = {0.0,0.0,0.0};

	while (it5.isNext())
	{
		auto key = it5.get();

		double tmp;

        tmp = 1.0f/g_vort.spacing(x)/2.0f*(g_vort.template get<vorticity>(key.move(x,1))[x] - g_vort.template get<vorticity>(key.move(x,-1))[x] ) +
			         	 	 1.0f/g_vort.spacing(y)/2.0f*(g_vort.template get<vorticity>(key.move(y,1))[y] - g_vort.template get<vorticity>(key.move(y,-1))[y] ) +
							 1.0f/g_vort.spacing(z)/2.0f*(g_vort.template get<vorticity>(key.move(z,1))[z] - g_vort.template get<vorticity>(key.move(z,-1))[z] );

        int_vort[x] += g_vort.template get<vorticity>(key)[x];
        int_vort[y] += g_vort.template get<vorticity>(key)[y];
        int_vort[z] += g_vort.template get<vorticity>(key)[z];

        if (tmp > max_vort)
        	max_vort = tmp;

		++it5;
	}

	Vcluster<> & v_cl = create_vcluster();
	v_cl.max(max_vort);
	v_cl.sum(int_vort[0]);
	v_cl.sum(int_vort[1]);
	v_cl.sum(int_vort[2]);
	v_cl.execute();

	if (v_cl.getProcessUnitID() == 0)
	{std::cout << "Max div for vorticity " << max_vort << "   Integral: " << int_vort[0] << "  " << int_vort[1] << "   " << int_vort[2] << std::endl;}

	//! \cond [sanity_int_div] \endcond
}

/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Step 1: Initialization of the vortex ring # {#vic_ring_init}
 *
 * In this function we initialize the vortex ring. The vortex ring is
 * initialized accordingly to these formula.
 *
 * \f$ w(t = 0) =  \frac{\Gamma}{\pi \sigma^{2}} e^{-(s/ \sigma)^2} \f$
 *
 * \f$ s^2 = (z-z_c)^{2} + ((x-x_c)^2 + (y-y_c)^2 - R^2) \f$
 *
 * \f$ \Gamma = \nu Re \f$
 *
 * With this initialization the vortex ring look like the one in figure
 *
 * \image html int_vortex_arrow_small.jpg "Vortex ring initialization the arrow indicate the direction where the vortex point while the colors indicate the magnitude from blue (low) to red (high)"
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp init_vort
 *
 *
 *
 */

//! \cond [init_vort] \endcond

/*
 * gr is the grid where we are initializing the vortex ring
 * domain is the simulation domain
 *
 */
void init_ring(grid_type & gr, const Box<3,float> & domain)
{
	// To add some noise to the vortex ring we create two random
	// vector
	constexpr int nk = 32;

	float ak[nk];
	float bk[nk];

	for (size_t i = 0 ; i < nk ; i++)
	{
	     ak[i] = rand()/RAND_MAX;
	     bk[i] = rand()/RAND_MAX;
	}

	// We calculate the circuitation gamma
	float gamma = nu * tgtre;
	float rinv2 = 1.0f/(sigma*sigma);
	float max_vorticity = gamma*rinv2/M_PI;

	// We go through the grid to initialize the vortex
	auto it = gr.getDomainIterator();

	while (it.isNext())
	{
		auto key_d = it.get();
		auto key = it.getGKey(key_d);

        float tx = (key.get(x)-2)*gr.spacing(x) + domain.getLow(x);
        float ty = (key.get(y)-2)*gr.spacing(y) + domain.getLow(y);
        float tz = (key.get(z)-2)*gr.spacing(z) + domain.getLow(z);
        float theta1 = atan2((ty-domain.getHigh(1)/2.0f),(tz-domain.getHigh(2)/2.0f));


        float rad1r  = sqrt((ty-domain.getHigh(1)/2.0f)*(ty-domain.getHigh(1)/2.0f) + (tz-domain.getHigh(2)/2.0f)*(tz-domain.getHigh(2)/2.0f)) - ringr1;
        float rad1t = tx - 1.0f;
        float rad1sq = rad1r*rad1r + rad1t*rad1t;
        float radstr = -exp(-rad1sq*rinv2)*rinv2*gamma/M_PI;
        gr.template get<vorticity>(key_d)[x] = 0.0f;
        gr.template get<vorticity>(key_d)[y] = -radstr * cos(theta1);
        gr.template get<vorticity>(key_d)[z] = radstr * sin(theta1);

        // kill the axis term

        float rad1r_  = sqrt((ty-domain.getHigh(1)/2.0f)*(ty-domain.getHigh(1)/2.0f) + (tz-domain.getHigh(2)/2.0f)*(tz-domain.getHigh(2)/2.0f)) + ringr1;
        float rad1sqTILDA = rad1sq*rinv2;
        radstr = exp(-rad1sq*rinv2)*rinv2*gamma/M_PI;
        gr.template get<vorticity>(key_d)[x] = 0.0f;
        gr.template get<vorticity>(key_d)[y] = -radstr * cos(theta1);
        gr.template get<vorticity>(key_d)[z] = radstr * sin(theta1);

		++it;
	}
}

//! \cond [init_vort] \endcond

//! \cond [poisson_syseq] \endcond

// Specification of the poisson equation for the helmotz-hodge projection
// 3D (dims = 3). The field is a scalar value (nvar = 1), bournary are periodic
// type of the the space is float. The grid type that store \psi
// The others indicate which kind of Matrix to use to construct the linear system and
// which kind of vector to construct for the right hand side. Here we use a PETSC Sparse Matrix
// and PETSC vector. NORMAL_GRID indicate that is a standard grid (non-staggered)
struct poisson_nn_helm
{
		//! 3D Stystem
        static const unsigned int dims = 3;
        //! We are solving for \psi that is a scalar
        static const unsigned int nvar = 1;
        //! Boundary conditions
        static const bool boundary[];
        //! type of the spatial coordinates
        typedef float stype;
        //! grid that store \psi
        typedef grid_dist_id<3,float,aggregate<float>> b_grid;
        //! Sparse matrix used to sove the linear system (we use PETSC)
        typedef SparseMatrix<double,int,PETSC_BASE> SparseMatrix_type;
        //! Vector to solve the system (PETSC)
        typedef Vector<double,PETSC_BASE> Vector_type;
        //! It is a normal grid
        static const int grid_type = NORMAL_GRID;
};

//! boundary conditions are PERIODIC
const bool poisson_nn_helm::boundary[] = {PERIODIC,PERIODIC,PERIODIC};

//! \cond [poisson_syseq] \endcond


/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Step 2: Helmotz-hodge projection # {#vic_hlm_proj}
 *
 * The Helnotz hodge projection is required in order to make the vorticity divergent
 * free. The Helmotz-holde projection work in this way. A field can be divided into
 * a curl-free part and a divergent-free part.
 *
 * \f$ w = w_{rot} + w_{div} \f$
 *
 * with
 *
 * \f$ \vec \nabla \times w_{rot} = 0 \f$
 *
 * \f$  \nabla \cdot w_{div} = 0 \f$
 *
 * To have a vorticity divergent free we have to get the component (3) \f$w_{div} = w - w_{rot}\f$.
 * In particular it hold
 *
 * \f$ \nabla \cdot w = \nabla \cdot w_{rot} \f$
 *
 * Bacause \f$ \vec \nabla \times w_{rot} = 0 \f$ we can introduce a field \f$ \psi \f$
 * such that
 *
 * (2) \f$ w_{rot} = \vec \nabla \psi \f$
 *
 * Doing the  \f$  \nabla \cdot \vec \nabla \psi \f$ we obtain
 *
 * \f$ \nabla \cdot \vec \nabla \psi = \nabla^{2} \psi = \nabla \cdot w_{rot} = \nabla \cdot w \f$
 *
 * so we lead to this equation
 *
 * (1) \f$ \nabla^{2} \psi = \nabla \cdot w  \f$
 *
 * Solving the equation for \f$ \psi \f$ we can obtain \f$ w_{rot} \f$ doing the gradient of \f$ \psi \f$
 * and finally correct \f$ w \f$ obtaining \f$ w_{div} \f$
 *
 * The **helmotz_hodge_projection** function do this correction to the vorticity
 *
 * In particular it solve the equation (1) it calculate \f$ w_{rot} \f$
 * using (2) and correct the vorticity using using (3)
 *
 *
 * ## Poisson equation ##
 *
 * To solve a poisson equation on a grid using finite-difference, we need to create
 * an object that carry information about the system of equations
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp poisson_syseq
 *
 * Once created this object we can define the equation we are trying to solve.
 * In particular the code below define the left-hand-side of the equation \f$ \nabla^{2} \psi \f$
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp poisson_obj_eq
 *
 * Before to construct the linear system we also calculate the divergence of the
 * vorticity \f$ \nabla \cdot w \f$ that will be the right-hand-side
 * of the equation
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp calc_div_vort
 *
 * Finally we can create the object FDScheme using the object **poisson_nn_helm**
 * as template variable. In addition to the constructor we have to specify the maximum extension of the stencil, the domain and the
 * grid that will store the result. At this point we can impose an equation to construct
 * our SparseMatrix. In this example we are imposing the poisson equation with right hand
 * side equal to the divergence of vorticity (note: to avoid to create another field we use
 *  \f$ \psi \f$ to preliminary store the divergence of the vorticity). Imposing the
 *  equations produce an invertible SparseMatrix **A** and a right-hand-side Vector **b**.
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp create_fdscheme
 *
 * Because we need \f$ x = A^{-1}b \f$. We have to invert and solve a linear system.
 * In this case we use the Conjugate-gradient-Method an iterative solver. Such method
 * is controlled by two parameters. One is the tollerance that determine when the
 * method is converged, the second one is the maximum number of iterations to avoid that
 * the method go into infinite loop. After we set the parameters of the solver we can the
 * the solution **x**. Finaly we copy back the solution **x** into the grid \f$ \psi \f$.
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp solve_petsc
 *
 * ### Note ###
 *
 * Because we are solving the poisson equation in periodic boundary conditions the Matrix has
 * determinat equal to zero. This mean that \f$ \psi \f$ has no unique solution (if it has one).
 * In order to recover one, we have to ensure that the integral of the righ hand side or vorticity
 * is zero. (In our case is the case). We have to ensure that across time the integral of the
 * vorticity is conserved. (In our case is the case if we consider the \f$ \nu = 0 \f$ and \f$
 * \nabla \cdot w = 0 \f$ we can rewrite (5) in a conservative way \f$  \frac{Dw}{dt} = div(w \otimes v) \f$ ).
 * Is also good to notice that the solution that you get is the one with \f$ \int w  = 0 \f$
 *
 *
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp solve_petsc
 *
 * ## Correction ## {#vort_correction}
 *
 * After we got our solution for \f$ \psi \f$ we can calculate the correction of the vorticity
 * doing the gradient of \f$ \psi \f$.
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp vort_correction
 *
 * We also do a sanity check and we control that the vorticity remain
 * divergent-free. Getting the maximum value of the divergence and printing out
 * its value
 *
 *
 */

/*
 * gr vorticity grid where we apply the correction
 * domain simulation domain
 *
 */
void helmotz_hodge_projection(grid_type & gr, const Box<3,float> & domain, petsc_solver<double> & solver, petsc_solver<double>::return_type & x_ , bool init)
{
	///////////////////////////////////////////////////////////////

	 //! \cond [poisson_obj_eq] \endcond

	// Convenient constant
	constexpr unsigned int phi = 0;

	// We define a field phi_f
	typedef Field<phi,poisson_nn_helm> phi_f;

	// We assemble the poisson equation doing the
	// poisson of the Laplacian of phi using Laplacian
	// central scheme (where the both derivative use central scheme, not
	// forward composed backward like usually)
	typedef Lap<phi_f,poisson_nn_helm,CENTRAL_SYM> poisson;

	//! \cond [poisson_obj_eq] \endcond

	//! \cond [calc_div_vort] \endcond

	// ghost get
	gr.template ghost_get<vorticity>();

	// ghost size of the psi function
    Ghost<3,long int> g(2);

	// Here we create a distributed grid to store the result of the helmotz projection
	grid_dist_id<3,float,aggregate<float>> psi(gr.getDecomposition(),gr.getGridInfo().getSize(),g);

	// Calculate the divergence of the vortex field
	auto it = gr.getDomainIterator();	

	while (it.isNext())
	{
		auto key = it.get();

        psi.template get<phi>(key) = 1.0f/gr.spacing(x)/2.0f*(gr.template get<vorticity>(key.move(x,1))[x] - gr.template get<vorticity>(key.move(x,-1))[x] ) +
			         	 	 1.0f/gr.spacing(y)/2.0f*(gr.template get<vorticity>(key.move(y,1))[y] - gr.template get<vorticity>(key.move(y,-1))[y] ) +
							 1.0f/gr.spacing(z)/2.0f*(gr.template get<vorticity>(key.move(z,1))[z] - gr.template get<vorticity>(key.move(z,-1))[z] );

		++it;
	}

	//! \cond [calc_div_vort] \endcond


	calc_and_print_max_div_and_int(gr);


	//! \cond [create_fdscheme] \endcond

	// In order to create a matrix that represent the poisson equation we have to indicate
	// we have to indicate the maximum extension of the stencil and we we need an extra layer
	// of points in case we have to impose particular boundary conditions.
	Ghost<3,long int> stencil_max(2);

	// Here we define our Finite difference disctetization scheme object
	FDScheme<poisson_nn_helm> fd(stencil_max, domain, psi);

	// impose the poisson equation, using right hand side psi on the full domain (psi.getDomainIterator)
	// the template paramereter instead define which property of phi carry the righ-hand-side
	// in this case phi has only one property, so the property 0 carry the right hand side
	fd.template impose_dit<0>(poisson(),psi,psi.getDomainIterator());

	//! \cond [create_fdscheme] \endcond

	//! \cond [solve_petsc] \endcond

	timer tm_solve;
	if (init == true)
	{
		// Set the conjugate-gradient as method to solve the linear solver
		solver.setSolver(KSPBCGS);

		// Set the absolute tolerance to determine that the method is converged
		solver.setAbsTol(0.0001);

		// Set the maximum number of iterations
		solver.setMaxIter(500);

#ifdef USE_MULTIGRID

        solver.setPreconditioner(PCHYPRE_BOOMERAMG);
        solver.setPreconditionerAMG_nl(6);
        solver.setPreconditionerAMG_maxit(1);
        solver.setPreconditionerAMG_relax("SOR/Jacobi");
        solver.setPreconditionerAMG_cycleType("V",0,4);
        solver.setPreconditionerAMG_coarsenNodalType(0);
        solver.setPreconditionerAMG_coarsen("HMIS");


#endif

		tm_solve.start();

		// Give to the solver A and b, return x, the solution
		solver.solve(fd.getA(),x_,fd.getB());

		tm_solve.stop();
	}
	else
	{
		tm_solve.start();
		solver.solve(x_,fd.getB());
		tm_solve.stop();
	}

	// copy the solution x to the grid psi
	fd.template copy<phi>(x_,psi);

	//! \cond [solve_petsc] \endcond

	//! \cond [vort_correction] \endcond

	psi.template ghost_get<phi>();

	// Correct the vorticity to make it divergence free

	auto it2 = gr.getDomainIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		gr.template get<vorticity>(key)[x] -= 1.0f/2.0f/psi.spacing(x)*(psi.template get<phi>(key.move(x,1))-psi.template get<phi>(key.move(x,-1)));
		gr.template get<vorticity>(key)[y] -= 1.0f/2.0f/psi.spacing(y)*(psi.template get<phi>(key.move(y,1))-psi.template get<phi>(key.move(y,-1)));
		gr.template get<vorticity>(key)[z] -= 1.0f/2.0f/psi.spacing(z)*(psi.template get<phi>(key.move(z,1))-psi.template get<phi>(key.move(z,-1)));

		++it2;
	}

	//! \cond [vort_correction] \endcond

	calc_and_print_max_div_and_int(gr);
}


/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Step 3: Remeshing vorticity # {#vic_remesh_vort}
 *
 * After that we initialized the vorticity on the grid, we initialize the particles
 * in a grid like position and we interpolate the vorticity on particles. Because
 * of the particles position being in a grid-like position and the symmetry of the
 * interpolation kernels, the re-mesh step simply reduce to initialize the particle
 * in a grid like position and assign the property vorticity of the particles equal to the
 * grid vorticity.
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp remesh_part
 *
 */

//! \cond [remesh_part] \endcond

void remesh(particles_type & vd, grid_type & gr,Box<3,float> & domain)
{
	// Remove all particles because we reinitialize in a grid like position
	vd.clear();

	// Get a grid iterator
	auto git = vd.getGridIterator(gr.getGridInfo().getSize());

	// For each grid point
	while (git.isNext())
	{
		// Get the grid point in global coordinates (key). p is in local coordinates
		auto p = git.get();
		auto key = git.get_dist();

		// Add the particle
		vd.add();

		// Assign the position to the particle in a grid like way
		vd.getLastPos()[x] = gr.spacing(x)*p.get(x) + domain.getLow(x);
		vd.getLastPos()[y] = gr.spacing(y)*p.get(y) + domain.getLow(y);
		vd.getLastPos()[z] = gr.spacing(z)*p.get(z) + domain.getLow(z);

		// Initialize the vorticity
		vd.template getLastProp<vorticity>()[x] = gr.template get<vorticity>(key)[x];
		vd.template getLastProp<vorticity>()[y] = gr.template get<vorticity>(key)[y];
		vd.template getLastProp<vorticity>()[z] = gr.template get<vorticity>(key)[z];

		// next grid point
		++git;
	}

	// redistribute the particles
	vd.map();
}

//! \cond [remesh_part] \endcond


/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Step 5: Compute velocity from vorticity # {#vic_vel_from_vort}
 *
 * Computing the velocity from vorticity is done in the following way. Given
 *
 * \f$ \vec \nabla \times u = -w \f$
 *
 * We intrododuce the stream line function defined as
 *
 * \f$ \nabla \times \phi = u \f$ (7)
 *
 * \f$ \nabla \cdot \phi = 0 \f$
 *
 * We obtain
 *
 * \f$ \nabla \times \nabla \times \phi = -w = \vec \nabla (\nabla \cdot  \phi) - \nabla^{2} \phi  \f$
 *
 * Because the divergence of \f$ \phi \f$ is defined to be zero we have
 *
 * \f$ \nabla^{2} \phi = w \f$
 *
 * The velocity can be recovered by the equation (7)
 *
 * Putting into code what explained before, we again generate a poisson
 * object
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp poisson_obj_eq
 *
 * In order to calculate the velocity out of the vorticity, we solve a poisson
 * equation like we did in helmotz-projection equation, but we do it for each
 * component \f$ i \f$ of the vorticity. Qnce we have the solution in **psi_s**
 * we copy the result back into the grid **gr_ps**. We than calculate the
 * quality of the solution printing the norm infinity of the residual and
 * finally we save in the grid vector vield **phi_v** the compinent \f$ i \f$
 * (Copy from phi_s to phi_v is necessary because in phi_s is not a grid
 * and cannot be used as a grid like object)
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp solve_poisson_comp
 *
 * We save the component \f$ i \f$ of \f$ \phi \f$ into **phi_v**
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp copy_to_phi_v
 *
 * Once we filled phi_v we can implement (7) and calculate the curl of **phi_v**
 * to recover the velocity v
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp curl_phi_v
 *
 *
 */

void comp_vel(Box<3,float> & domain, grid_type & g_vort,grid_type & g_vel, petsc_solver<double>::return_type (& phi_s)[3],petsc_solver<double> & solver)
{
	//! \cond [poisson_obj_eq] \endcond

	// Convenient constant
	constexpr unsigned int phi = 0;

	// We define a field phi_f
	typedef Field<phi,poisson_nn_helm> phi_f;

	// We assemble the poisson equation doing the
	// poisson of the Laplacian of phi using Laplacian
	// central scheme (where the both derivative use central scheme, not
	// forward composed backward like usually)
	typedef Lap<phi_f,poisson_nn_helm,CENTRAL_SYM> poisson;

	//! \cond [poisson_obj_eq] \endcond

	// Maximum size of the stencil
	Ghost<3,long int> g(2);

	// Here we create a distributed grid to store the result
	grid_dist_id<3,float,aggregate<float>> gr_ps(g_vort.getDecomposition(),g_vort.getGridInfo().getSize(),g);
	grid_dist_id<3,float,aggregate<float[3]>> phi_v(g_vort.getDecomposition(),g_vort.getGridInfo().getSize(),g);

	// We calculate and print the maximum divergence of the vorticity
	// and the integral of it
	calc_and_print_max_div_and_int(g_vort);

	// For each component solve a poisson
	for (size_t i = 0 ; i < 3 ; i++)
	{
		//! \cond [solve_poisson_comp] \endcond

		// Copy the vorticity component i in gr_ps
		auto it2 = gr_ps.getDomainIterator();

		// calculate the velocity from the curl of phi
		while (it2.isNext())
		{
			auto key = it2.get();

			// copy
			gr_ps.get<vorticity>(key) = g_vort.template get<vorticity>(key)[i];

			++it2;
		}

		// Maximum size of the stencil
		Ghost<3,long int> stencil_max(2);

		// Finite difference scheme
		FDScheme<poisson_nn_helm> fd(stencil_max, domain, gr_ps);

		// impose the poisson equation using gr_ps = vorticity for the right-hand-side (on the full domain)
		fd.template impose_dit<phi>(poisson(),gr_ps,gr_ps.getDomainIterator());

		solver.setAbsTol(0.01);

		// Get the vector that represent the right-hand-side
		Vector<double,PETSC_BASE> & b = fd.getB();

		timer tm_solve;
		tm_solve.start();

		// Give to the solver A and b in this case we are giving
		// an intitial guess phi_s calculated in the previous
		// time step
		solver.solve(phi_s[i],b);

		tm_solve.stop();

		// Calculate the residual

		solError serr;
		serr = solver.get_residual_error(phi_s[i],b);

		Vcluster<> & v_cl = create_vcluster();
		if (v_cl.getProcessUnitID() == 0)
		{std::cout << "Solved component " << i << "  Error: " << serr.err_inf << std::endl;}

		// copy the solution to grid
		fd.template copy<phi>(phi_s[i],gr_ps);

		//! \cond [solve_poisson_comp] \endcond

		//! \cond [copy_to_phi_v] \endcond

		auto it3 = gr_ps.getDomainIterator();

		// calculate the velocity from the curl of phi
		while (it3.isNext())
		{
			auto key = it3.get();

			phi_v.get<velocity>(key)[i] = gr_ps.get<phi>(key);

			++it3;
		}

		//! \cond [copy_to_phi_v] \endcond
	}

	//! \cond [curl_phi_v] \endcond

	phi_v.ghost_get<phi>();

	float inv_dx = 0.5f / g_vort.spacing(0);
	float inv_dy = 0.5f / g_vort.spacing(1);
	float inv_dz = 0.5f / g_vort.spacing(2);

	auto it3 = phi_v.getDomainIterator();

	// calculate the velocity from the curl of phi
	while (it3.isNext())
	{
		auto key = it3.get();

		float phi_xy = (phi_v.get<phi>(key.move(y,1))[x] - phi_v.get<phi>(key.move(y,-1))[x])*inv_dy;
		float phi_xz = (phi_v.get<phi>(key.move(z,1))[x] - phi_v.get<phi>(key.move(z,-1))[x])*inv_dz;
		float phi_yx = (phi_v.get<phi>(key.move(x,1))[y] - phi_v.get<phi>(key.move(x,-1))[y])*inv_dx;
		float phi_yz = (phi_v.get<phi>(key.move(z,1))[y] - phi_v.get<phi>(key.move(z,-1))[y])*inv_dz;
		float phi_zx = (phi_v.get<phi>(key.move(x,1))[z] - phi_v.get<phi>(key.move(x,-1))[z])*inv_dx;
		float phi_zy = (phi_v.get<phi>(key.move(y,1))[z] - phi_v.get<phi>(key.move(y,-1))[z])*inv_dy;

		g_vel.template get<velocity>(key)[x] = phi_zy - phi_yz;
		g_vel.template get<velocity>(key)[y] = phi_xz - phi_zx;
		g_vel.template get<velocity>(key)[z] = phi_yx - phi_xy;

		++it3;
	}

	//! \cond [curl_phi_v] \endcond
}


template<unsigned int prp> void set_zero(grid_type & gr)
{
	auto it = gr.getDomainGhostIterator();

	// calculate the velocity from the curl of phi
	while (it.isNext())
	{
		auto key = it.get();

		gr.template get<prp>(key)[0] = 0.0;
		gr.template get<prp>(key)[1] = 0.0;
		gr.template get<prp>(key)[2] = 0.0;

		++it;
	}
}


template<unsigned int prp> void set_zero(particles_type & vd)
{
	auto it = vd.getDomainIterator();

	// calculate the velocity from the curl of phi
	while (it.isNext())
	{
		auto key = it.get();

		vd.template getProp<prp>(key)[0] = 0.0;
		vd.template getProp<prp>(key)[1] = 0.0;
		vd.template getProp<prp>(key)[2] = 0.0;

		++it;
	}
}

/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Step 6: Compute right hand side # {#vic_rhs_calc}
 *
 * Computing the right hand side is performed calculating the term
 * \f$ (w \cdot \nabla) u \f$. For the nabla operator we use second
 * order finite difference central scheme. The commented part is the
 * term \f$ \nu \nabla^{2} w \f$ that we said to neglect
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp calc_rhs
 *
 */

//! \cond [calc_rhs] \endcond

// Calculate the right hand side of the vorticity formulation
template<typename grid> void calc_rhs(grid & g_vort, grid & g_vel, grid & g_dwp)
{
	// usefull constant
	constexpr int rhs = 0;

	// calculate several pre-factors for the stencil finite
	// difference
	float fac1 = 2.0f*nu/(g_vort.spacing(0)*g_vort.spacing(0));
	float fac2 = 2.0f*nu/(g_vort.spacing(1)*g_vort.spacing(1));
	float fac3 = 2.0f*nu/(g_vort.spacing(2)*g_vort.spacing(2));

	float fac4 = 0.5f/(g_vort.spacing(0));
	float fac5 = 0.5f/(g_vort.spacing(1));
	float fac6 = 0.5f/(g_vort.spacing(2));

	auto it = g_dwp.getDomainIterator();

	g_vort.template ghost_get<vorticity>();

	while (it.isNext())
	{
		auto key = it.get();

		g_dwp.template get<rhs>(key)[x] = fac1*(g_vort.template get<vorticity>(key.move(x,1))[x]+g_vort.template get<vorticity>(key.move(x,-1))[x])+
				                          fac2*(g_vort.template get<vorticity>(key.move(y,1))[x]+g_vort.template get<vorticity>(key.move(y,-1))[x])+
										  fac3*(g_vort.template get<vorticity>(key.move(z,1))[x]+g_vort.template get<vorticity>(key.move(z,-1))[x])-
										  2.0f*(fac1+fac2+fac3)*g_vort.template get<vorticity>(key)[x]+
										  fac4*g_vort.template get<vorticity>(key)[x]*
										  (g_vel.template get<velocity>(key.move(x,1))[x] - g_vel.template get<velocity>(key.move(x,-1))[x])+
										  fac5*g_vort.template get<vorticity>(key)[y]*
										  (g_vel.template get<velocity>(key.move(y,1))[x] - g_vel.template get<velocity>(key.move(y,-1))[x])+
										  fac6*g_vort.template get<vorticity>(key)[z]*
										  (g_vel.template get<velocity>(key.move(z,1))[x] - g_vel.template get<velocity>(key.move(z,-1))[x]);

		g_dwp.template get<rhs>(key)[y] = fac1*(g_vort.template get<vorticity>(key.move(x,1))[y]+g_vort.template get<vorticity>(key.move(x,-1))[y])+
				                          fac2*(g_vort.template get<vorticity>(key.move(y,1))[y]+g_vort.template get<vorticity>(key.move(y,-1))[y])+
										  fac3*(g_vort.template get<vorticity>(key.move(z,1))[y]+g_vort.template get<vorticity>(key.move(z,-1))[y])-
										  2.0f*(fac1+fac2+fac3)*g_vort.template get<vorticity>(key)[y]+
										  fac4*g_vort.template get<vorticity>(key)[x]*
										  (g_vel.template get<velocity>(key.move(x,1))[y] - g_vel.template get<velocity>(key.move(x,-1))[y])+
										  fac5*g_vort.template get<vorticity>(key)[y]*
										  (g_vel.template get<velocity>(key.move(y,1))[y] - g_vel.template get<velocity>(key.move(y,-1))[y])+
										  fac6*g_vort.template get<vorticity>(key)[z]*
										  (g_vel.template get<velocity>(key.move(z,1))[y] - g_vel.template get<velocity>(key.move(z,-1))[y]);


		g_dwp.template get<rhs>(key)[z] = fac1*(g_vort.template get<vorticity>(key.move(x,1))[z]+g_vort.template get<vorticity>(key.move(x,-1))[z])+
				                          fac2*(g_vort.template get<vorticity>(key.move(y,1))[z]+g_vort.template get<vorticity>(key.move(y,-1))[z])+
										  fac3*(g_vort.template get<vorticity>(key.move(z,1))[z]+g_vort.template get<vorticity>(key.move(z,-1))[z])-
										  2.0f*(fac1+fac2+fac3)*g_vort.template get<vorticity>(key)[z]+
										  fac4*g_vort.template get<vorticity>(key)[x]*
										  (g_vel.template get<velocity>(key.move(x,1))[z] - g_vel.template get<velocity>(key.move(x,-1))[z])+
										  fac5*g_vort.template get<vorticity>(key)[y]*
										  (g_vel.template get<velocity>(key.move(y,1))[z] - g_vel.template get<velocity>(key.move(y,-1))[z])+
										  fac6*g_vort.template get<vorticity>(key)[z]*
										  (g_vel.template get<velocity>(key.move(z,1))[z] - g_vel.template get<velocity>(key.move(z,-1))[z]);

		++it;
	}
}

//! \cond [calc_rhs] \endcond

/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Step 8: Runge-Kutta # {#vic_runge_kutta1}
 *
 * Here we do the first step of the runge kutta update. In particular we
 * update the vorticity and position of the particles. The right-hand-side
 * of the vorticity update is calculated on the grid and interpolated
 *  on the particles. The Runge-Kutta of order two
 * require the following update for the vorticity and position as first step
 *
 * \f$ \boldsymbol w = \boldsymbol w + \frac{1}{2} \boldsymbol {rhs} \delta t \f$
 *
 * \f$ \boldsymbol x = \boldsymbol x + \frac{1}{2} \boldsymbol u \delta t \f$
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp runge_kutta_1
 *
 */

//! \cond [runge_kutta_1] \endcond

void rk_step1(particles_type & particles)
{
	auto it = particles.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		// Save the old vorticity
		particles.getProp<old_vort>(key)[x] = particles.getProp<vorticity>(key)[x];
		particles.getProp<old_vort>(key)[y] = particles.getProp<vorticity>(key)[y];
		particles.getProp<old_vort>(key)[z] = particles.getProp<vorticity>(key)[z];

		particles.getProp<vorticity>(key)[x] = particles.getProp<vorticity>(key)[x] + 0.5f * dt * particles.getProp<rhs_part>(key)[x];
		particles.getProp<vorticity>(key)[y] = particles.getProp<vorticity>(key)[y] + 0.5f * dt * particles.getProp<rhs_part>(key)[y];
		particles.getProp<vorticity>(key)[z] = particles.getProp<vorticity>(key)[z] + 0.5f * dt * particles.getProp<rhs_part>(key)[z];

		++it;
	}

	auto it2 = particles.getDomainIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		// Save the old position
		particles.getProp<old_pos>(key)[x] = particles.getPos(key)[x];
		particles.getProp<old_pos>(key)[y] = particles.getPos(key)[y];
		particles.getProp<old_pos>(key)[z] = particles.getPos(key)[z];

		particles.getPos(key)[x] = particles.getPos(key)[x] + 0.5f * dt * particles.getProp<p_vel>(key)[x];
		particles.getPos(key)[y] = particles.getPos(key)[y] + 0.5f * dt * particles.getProp<p_vel>(key)[y];
		particles.getPos(key)[z] = particles.getPos(key)[z] + 0.5f * dt * particles.getProp<p_vel>(key)[z];

		++it2;
	}

	particles.map();
}

//! \cond [runge_kutta_1] \endcond

/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Step 13: Runge-Kutta # {#vic_runge_kutta2}
 *
 * Here we do the second step of the Runge-Kutta update. In particular we
 * update the vorticity and position of the particles. The right-hand-side
 * of the vorticity update is calculated on the grid and interpolated
 *  on the particles. The Runge-Kutta of order two
 * require the following update for the vorticity and position as first step
 *
 * \f$ \boldsymbol w = \boldsymbol w + \frac{1}{2} \boldsymbol {rhs} \delta t \f$
 *
 * \f$ \boldsymbol x = \boldsymbol x + \frac{1}{2} \boldsymbol u \delta t \f$
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp runge_kutta_2
 *
 */

//! \cond [runge_kutta_2] \endcond

void rk_step2(particles_type & particles)
{
	auto it = particles.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		particles.getProp<vorticity>(key)[x] = particles.getProp<old_vort>(key)[x] + dt * particles.getProp<rhs_part>(key)[x];
		particles.getProp<vorticity>(key)[y] = particles.getProp<old_vort>(key)[y] + dt * particles.getProp<rhs_part>(key)[y];
		particles.getProp<vorticity>(key)[z] = particles.getProp<old_vort>(key)[z] + dt * particles.getProp<rhs_part>(key)[z];

		++it;
	}

	auto it2 = particles.getDomainIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		particles.getPos(key)[x] = particles.getProp<old_pos>(key)[x] + dt * particles.getProp<p_vel>(key)[x];
		particles.getPos(key)[y] = particles.getProp<old_pos>(key)[y] + dt * particles.getProp<p_vel>(key)[y];
		particles.getPos(key)[z] = particles.getProp<old_pos>(key)[z] + dt * particles.getProp<p_vel>(key)[z];

		++it2;
	}

	particles.map();
}



//! \cond [runge_kutta_2] \endcond

/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Step 4-5-6-7: Do step # {#vic_do_step}
 *
 * The do step function assemble multiple steps some of them already explained.
 * First we interpolate the vorticity from particles to mesh
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp do_step_p2m
 *
 * than we calculate velocity out of vorticity and the right-hand-side
 * recalling step 5 and 6
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp step_56
 *
 * Finally we interpolate velocity and right-hand-side back to particles
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp inte_m2p
 *
 */

template<typename grid, typename vector> void do_step(vector & particles,
		                                              grid & g_vort,
													  grid & g_vel,
													  grid & g_dvort,
													  Box<3,float> & domain,
													  interpolate<particles_type,grid_type,mp4_kernel<float>> & inte,
													  petsc_solver<double>::return_type (& phi_s)[3],
													  petsc_solver<double> & solver)
{
	constexpr int rhs = 0;

	//! \cond [do_step_p2m] \endcond

	set_zero<vorticity>(g_vort);
	inte.template p2m<vorticity,vorticity>(particles,g_vort);

	g_vort.template ghost_put<add_,vorticity>();

	//! \cond [do_step_p2m] \endcond

	//! \cond [step_56] \endcond

	// Calculate velocity from vorticity
	comp_vel(domain,g_vort,g_vel,phi_s,solver);
	calc_rhs(g_vort,g_vel,g_dvort);

	//! \cond [step_56] \endcond

	//! \cond [inte_m2p] \endcond

	g_dvort.template ghost_get<rhs>();
	g_vel.template ghost_get<velocity>();
	set_zero<rhs_part>(particles);
	set_zero<p_vel>(particles);
	inte.template m2p<rhs,rhs_part>(g_dvort,particles);
	inte.template m2p<velocity,p_vel>(g_vel,particles);

	//! \cond [inte_m2p] \endcond
}

template<typename vector, typename grid> void check_point_and_save(vector & particles,
																   grid & g_vort,
																   grid & g_vel,
																   grid & g_dvort,
																   size_t i)
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() < 24)
	{
		particles.write("part_out_" + std::to_string(i),VTK_WRITER | FORMAT_BINARY);
		g_vort.template ghost_get<vorticity>();
		g_vel.template ghost_get<velocity>();
		g_vel.write_frame("grid_velocity",i, VTK_WRITER | FORMAT_BINARY);
		g_vort.write_frame("grid_vorticity",i, VTK_WRITER | FORMAT_BINARY);
	}

	// In order to reduce the size of the saved data we apply a threshold.
	// We only save particles with vorticity higher than 0.1
	vector_dist<3,float,aggregate<float>> part_save(particles.getDecomposition(),0);

	auto it_s = particles.getDomainIterator();

	while (it_s.isNext())
	{
		auto p = it_s.get();

		float vort_magn = sqrt(particles.template getProp<vorticity>(p)[0] * particles.template getProp<vorticity>(p)[0] +
							   particles.template getProp<vorticity>(p)[1] * particles.template getProp<vorticity>(p)[1] +
							   particles.template getProp<vorticity>(p)[2] * particles.template getProp<vorticity>(p)[2]);

		if (vort_magn < 0.1)
		{
			++it_s;
			continue;
		}

		part_save.add();

		part_save.getLastPos()[0] = particles.getPos(p)[0];
		part_save.getLastPos()[1] = particles.getPos(p)[1];
		part_save.getLastPos()[2] = particles.getPos(p)[2];

		part_save.template getLastProp<0>() = vort_magn;

		++it_s;
	}

	part_save.map();

	// Save and HDF5 file for checkpoint restart
	particles.save("check_point");
}

/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Main # {#vic_main}
 *
 * The main function as usual call the function **openfpm_init** it define
 * the domain where the simulation take place. The ghost size of the grid
 * the size of the grid in grid units on each direction, the periodicity
 * of the domain, in this case PERIODIC in each dimension and we create
 * our basic data structure for the simulation. A grid for the vorticity
 * **g_vort** a grid for the velocity **g_vel** a grid for the right-hand-side
 * of the vorticity update, and the **particles** vector. Additionally
 * we define the data structure **phi_s[3]** that store the velocity solution
 * in the previous time-step. keeping track of the previous solution for
 * the velocity help the interative-solver to find the solution more quickly.
 * Using the old velocity configuration as initial guess the solver will
 * converge in few iterations refining the old one.
 *
 */
int main(int argc, char* argv[])
{
	// Initialize
	openfpm_init(&argc,&argv);
	{
	// Domain, a rectangle
	// For the grid 1600x400x400 use
	// Box<3,float> domain({0.0,0.0,0.0},{22.0,5.57,5.57});
	Box<3,float> domain({0.0,0.0,0.0},{3.57,3.57,3.57});

	// Ghost (Not important in this case but required)
	Ghost<3,long int> g(2);

	// Grid points on x=128 y=64 z=64
	// if we use Re = 7500
	//  long int sz[] = {1600,400,400};
	long int sz[] = {96,96,96};
	size_t szu[] = {(size_t)sz[0],(size_t)sz[1],(size_t)sz[2]};

	periodicity<3> bc = {{PERIODIC,PERIODIC,PERIODIC}};

	grid_type g_vort(szu,domain,g,bc);
	grid_type g_vel(g_vort.getDecomposition(),szu,g);
	grid_type g_dvort(g_vort.getDecomposition(),szu,g);
	particles_type particles(g_vort.getDecomposition(),0);

	// It store the solution to compute velocity
	// It is used as initial guess every time we call the solver
	Vector<double,PETSC_BASE> phi_s[3];
	Vector<double,PETSC_BASE> x_;

	// Parallel object
	Vcluster<> & v_cl = create_vcluster();

	// print out the spacing of the grid
	if (v_cl.getProcessUnitID() == 0)
	{std::cout << "SPACING " << g_vort.spacing(0) << " " << g_vort.spacing(1) << " " << g_vort.spacing(2) << std::endl;}

	// initialize the ring step 1
	init_ring(g_vort,domain);

	x_.resize(g_vort.size(),g_vort.getLocalDomainSize());
	x_.setZero();

	// Create a PETSC solver to get the solution x
	petsc_solver<double> solver;

	// Do the helmotz projection step 2
	helmotz_hodge_projection(g_vort,domain,solver,x_,true);

	// initialize the particle on a mesh step 3
	remesh(particles,g_vort,domain);

	// Get the total number of particles
	size_t tot_part = particles.size_local();
	v_cl.sum(tot_part);
	v_cl.execute();

	// Now we set the size of phi_s
	phi_s[0].resize(g_vort.size(),g_vort.getLocalDomainSize());
	phi_s[1].resize(g_vort.size(),g_vort.getLocalDomainSize());
	phi_s[2].resize(g_vort.size(),g_vort.getLocalDomainSize());

	// And we initially set it to zero
	phi_s[0].setZero();
	phi_s[1].setZero();
	phi_s[2].setZero();

	// We create the interpolation object outside the loop cycles
	// to avoid to recreate it at every cycle. Internally this object
	// create some data-structure to optimize the conversion particle
	// position to sub-domain. If there are no re-balancing it is safe
	// to reuse-it
	interpolate<particles_type,grid_type,mp4_kernel<float>> inte(particles,g_vort);

	// With more than 24 core we rely on the HDF5 checkpoint restart file
	if (v_cl.getProcessingUnits() < 24)
	{g_vort.write("grid_vorticity_init", VTK_WRITER | FORMAT_BINARY);}


	// Before entering the loop we check if we want to restart from
	// a previous simulation
	size_t i = 0;

	// if we have an argument
	if (argc > 1)
	{
		// take that argument
		std::string restart(argv[1]);

		// convert it into number
		i = std::stoi(restart);

		// load the file
		particles.load("check_point_" + std::to_string(i));

		// Print to inform that we are restarting from a
		// particular position
		if (v_cl.getProcessUnitID() == 0)
		{std::cout << "Restarting from " << i << std::endl;}
	}

	// Time Integration
	for ( ; i < nsteps ; i++)
	{
		// do step 4-5-6-7
		do_step(particles,g_vort,g_vel,g_dvort,domain,inte,phi_s,solver);

		// do step 8
		rk_step1(particles);

		// do step 9-10-11-12
		do_step(particles,g_vort,g_vel,g_dvort,domain,inte,phi_s,solver);

		// do step 13
		rk_step2(particles);

		// so step 14
		set_zero<vorticity>(g_vort);
		inte.template p2m<vorticity,vorticity>(particles,g_vort);
		g_vort.template ghost_put<add_,vorticity>();

		// helmotz-hodge projection
		helmotz_hodge_projection(g_vort,domain,solver,x_,false);

		remesh(particles,g_vort,domain);

		// print the step number
		if (v_cl.getProcessUnitID() == 0)
		{std::cout << "Step " << i << std::endl;}

		// every 100 steps write the output
		if (i % 100 == 0)		{check_point_and_save(particles,g_vort,g_vel,g_dvort,i);}

	}
	}

	openfpm_finalize();
}

#else

int main(int argc, char* argv[])
{
        return 0;
}

#endif
