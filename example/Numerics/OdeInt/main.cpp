//
// Created by Abhinav Singh @absingh on 26/04/2021
//
/**
 * @file Odeint/main.cpp
 * @page Odeint Odeint
 *
 *  * We can use the Odeint library from boost with OpenFPM in two ways. The first example will cover a simple case where we control the time stepping ourselves (time loop of the simulation wriiten by us).
 * The 2nd Example will cover a way of using odeint where the time is forwarded by odeint along with an observer.
 *
 *  * In the example links below, we solve the follwing advection-diffusion PDE for two chemicals with concentration \f$\vec{C}(x,y,t)=(C_x(x,y,t),C_y(x,y,t))\f$ with a fixed velocity \f$\vec{V}(x,y)=(-ye^{-10(x^2+y^2)},xe^{-10(x^2+y^2)})\f$:
 *
 * @f[ \frac{\partial\vec{C}}{dt}=\vec{V}.\nabla_{\{x,y\}} \vec{C} + 0.1\Delta_{\{x,y\}} \vec{C} @f]
 * in 2d domain [-1,-1]*[1,1] with the boundary conditions on the walls as no-slip for velocity \f$\vec{V}=0\f$ and sink for the chemicals \f$\vec{C}=0\f$ for all time \f$t\f$.
 * Further, we start with the initial condition for the concentration as
 *   @f[\vec{C}=\begin{cases}
 *   (1,0)\text{ for } x=0,-0.5<y<0\\
 *   (0,1)\text{ for } x=0, 0<y<0.5\\
 *   (0,0)  \text{ for the rest of the domain}\\
 * \end{cases} @f]
 *
 *
 * @subpage Odeint_single_step
 * @subpage Odeint_multiple_steps
 *
 *
 */


/**
 * @page Odeint_single_step Step by step time integration with Odeint
 *

 * In this example, we perform time integration in a 2d domain of particles of a following partial differential equation:
 *
 * @f[ \frac{\partial\vec{C}}{dt}=\vec{V}.\nabla C + 0.1\Delta C @f]
 * in 2d domain [-1,-1]*[1,1] with a fixed velocity \f$\vec{V}(x,y)=(-ye^{-10(x^2+y^2)},xe^{-10(x^2+y^2)})\f$, and the boundary conditions on the walls as no-slip for velocity \f$\vec{V}=0\f$ and sink for the chemicals \f$\vec{C}=0\f$ for all time \f$t\f$.
 *  * Further, we start with the initial condition for the concentration as
 *   @f[\vec{C}=\begin{cases}
 *   (1,0)\text{ for } x=0,-0.5<y<0\\
 *   (0,1)\text{ for } x=0, 0<y<0.5\\
 *   (0,0)  \text{ for the rest of the domain}\\
 * \end{cases} @f]
 *
 * We do that by employing a Lagrangian frame of reference. Hence the PDE is transformed to a system of PDEs:
 *
 * @f[\begin{align}
 *   \frac{\partial\vec{X}}{dt}=\vec{V}\\
 *   \frac{\partial\vec{C}}{dt}=0.1\Delta_{\{x,y\}} \vec{C}
 * \end{align} @f]
 * Where \f$\vec{X}\f$ is the position vector. Hence, we have decoupled the moving of the particles from the evolution of the chemicals. As it can be expensive to recompute derivatives at every stage of a time integrator in a single step, we will integrate the PDEs with seperate techniques (Euler step for moving the particles).
 *
 * Output:
 * Time series data of the PDE Solution.
 *
 */

/**
 * @page Odeint_single_step Step by step time integration with Odeint
 *
 * ## Including the headers ## {#ode_c1_include}
 *
 * These are the header files that we need to include:
 *
 * @snippet example/Numerics/OdeInt/main.cpp Ode1Include
 *
 */
//! @cond [Ode1Include] @endcond
// Include Vector Expression,Vector Expressions for Subset,DCPSE,Odeint header files
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "DCPSE/DCPSE_op/DCPSE_op.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"
//! @cond [Ode1Include] @endcond


/**
 * @page Odeint_single_step Step by step time integration with Odeint
 *
 * ## Initialization of the global parameters## {#ode_c_init}
 *
 * We start with
 * * Initializing certain global parameteres we will use: such as x,y to refer to the dimensions 0 and 1. (Makes it easier to read equations in code)
 * * Creating empty pointers for coupling openfpm distributed vector with odeint. One for the entire distributed vector and another for the subset or bulk.
 * We seperate bulk and the entire distribution as it makes it easier to impose boundary conditions. (Which will be more apparant in ComputeRHS of the PDE)
 * Note that a subset expression always comes at the left hand side of a computation. (The semantics of the expressions is by denoting what we want to update from regular expressions)
 *
 * Creating aliases of the types of the datasructures we are going to use in OpenFPM.
 * Property_type as the type of properties we wish to use.
 * dist_vector_type as the 2d openfpm distributed vector type
 * dist_vector_type as the 2d openfpm distributed subset vector type
 *
 * @snippet example/Numerics/OdeInt/main.cpp Init1pm
 *
 */

//! @cond [Init1pm] @endcond

constexpr int x = 0;
constexpr int y = 1;

double dt;

void *PointerDistGlobal, *PointerDistSubset;

typedef aggregate<VectorS<2, double>, VectorS<2, double>, VectorS<2, double>> Property_type;
typedef vector_dist_ws<2, double, Property_type> dist_vector_type;
typedef vector_dist_subset<2, double, Property_type> dist_vector_subset_type;
//! @cond [Init1pm] @endcond

/**
 * @page Odeint_single_step Step by step time integration with Odeint
 *
 * ## Creating the RHS Functor## {#ode_c1_rhs}
 *
 * Odeint works with certain specific state_types.
 * We offer certain state types such as 'state_type_2d_ofp' for making openfpm work with odeint.
 *
 *
 * Now we create the RHS functor. Please refer to ODE_int for more detials.
 * Note that we have templated it with two operator types DXX and DYY as we need to compute Laplacian at each stage. We will pass the DCPSE operators to an instance of this functor.
 *
 * All RHS computations  needs to happen in the operator ().
 * Odeint expects the arguments here to be an input state_type X, an output state_tyoe dxdt and time t.
 * We pass on the openfpm distributed state types as
 * void operator()( const state_type_2d_ofp &X , state_type_2d_ofp &dxdt , const double t ) const
 *
 * Since we would like to use openfpm here. We cast back the global pointers created before to access the Openfpm distributed vector here.
 * (Note that these pointers needs to initialized in the main(). Further, 'state_type_2d_ofp' is a temporal structure, which means it does not have the ghost.
 * Hence we copy the current state back to the openfpm vector from the openfpm state type X.
 * We do our computations as required.
 * Then we copy back the output into the state_type dxdt)
 *
 *
 *
 * @snippet example/Numerics/OdeInt/main.cpp RHS1Functor
 *
 */
//! @cond [RHS1Functor] @endcond
template<typename DXX,typename DYY>
struct RHSFunctor
{
    //Intializing the operators
    DXX &Dxx;
    DYY &Dyy;
    //Constructor
    RHSFunctor(DXX &Dxx,DYY &Dyy):Dxx(Dxx),Dyy(Dyy)
    {}

    void operator()( const state_type_2d_ofp &X , state_type_2d_ofp &dxdt , const double t ) const
    {
        //Casting the pointers to OpenFPM vector distributions
        dist_vector_type &Particles= *(dist_vector_type *) PointerDistGlobal;
        dist_vector_subset_type &Particles_bulk= *(dist_vector_subset_type *) PointerDistSubset;

        //Aliasing the properties.
        auto C = getV<0>(Particles);
        auto C_bulk = getV<0>(Particles_bulk);
        auto dC = getV<2>(Particles);
        auto dC_bulk = getV<2>(Particles_bulk);

        //Since the state vector does not have enough information to compute derivates, we copy them back to the Particle Distribution.
        //In this way we can compute all the things that are required to be computed in each stage of the ode integrator.
        // Note that we are going to use bulk expressions as needed. Whenever we don't want to update the boundaries (Due to the boundary conditions).

        //These expressions only update the bulk values of C.
        C_bulk[x]=X.data.get<0>();
        C_bulk[y]=X.data.get<1>();
        Particles.ghost_get<0>();

        // We do the RHS computations for the Laplacian (Updating bulk only).
        dC_bulk[x] = 0.1*(Dxx(C[x])+Dyy(C[x]));
        dC_bulk[y] = 0.1*(Dxx(C[y])+Dyy(C[y]));

        //We copy back to the dxdt state_type for Odeint
        dxdt.data.get<0>()=dC[x];
        dxdt.data.get<1>()=dC[y];
    }
};
//! @cond [RHS1Functor] @endcond



/**
 * @page Odeint_single_step Step by step time integration with Odeint
 *
 * ## Initializating OpenFPM ## {#odeint_c_initmain}
 *
 * We start with
 * * Initializing OpenFPM
 *
 * @snippet example/Numerics/OdeInt/main.cpp init1Particles
 *
 */
//! @cond [init1Particles] @endcond
int main(int argc, char* argv[]) {
    //	initialize library
    openfpm_init(&argc, &argv);
    //! @cond [init1Particles] @endcond

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * @page Odeint_single_step Step by step time integration with Odeint
     *
     * ## Creating Particles and assigning subsets ## {#odeint_c_indices}
     *
     * We create a particle distribution we certain rCut for the domain [-1,-1] to [1,1].
     *
     * Also, we fill the initial concentration as C_1(x=0,y>0 & y<0.5,t=0)=1,C_2(x=0,y<0 & y>-0.5,t=0)=1 and 0 everywhere else.
     * @snippet example/Numerics/OdeInt/main.cpp init1Subset
     *
     */
    //! @cond [init1Subset] @endcond
    const size_t sz[2] = {41, 41};
    Box<2, double> box({-1, -1}, {1.0, 1.0});
    size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
    double spacing[2];
    spacing[0] = 2.0 / (sz[0] - 1);
    spacing[1] = 2.0 / (sz[1] - 1);
    double rCut = 3.1 * spacing[0];
    Ghost<2, double> ghost(rCut);

    dist_vector_type Particles(0, box, bc, ghost);
    Particles.setPropNames({"Concentration", "Velocity", "TempConcentration"});

    auto it = Particles.getGridIterator(sz);
    while (it.isNext()) {
        Particles.add();
        auto key = it.get();
        double x = -1.0 + key.get(0) * spacing[0];
        Particles.getLastPos()[0] = x;
        double y = -1.0 + key.get(1) * spacing[1];
        Particles.getLastPos()[1] = y;
        //We detect Boundaries. By labelling the subset to be 0 for bulk and 1 for the boundary. Subsets will be constructed later based on the label.
        if (x != -1.0 && x != 1.0 && y != -1.0 && y != 1) {
            Particles.getLastSubset(0);
        } else {
            Particles.getLastSubset(1);
        }
        // Here fill the Initial value of the concentration.
        if (x == 0.0 && y > -0.5 && y < 0) {
            Particles.template getLastProp<0>()[0] = 1.0;
            Particles.template getLastProp<0>()[1] = 0.0;
        } else if (x == 0.0 && y > 0 && y < 0.5) {
            Particles.template getLastProp<0>()[0] = 0.0;
            Particles.template getLastProp<0>()[1] = 1.0;
        } else {
            Particles.template getLastProp<0>()[0] = 0.0;
            Particles.template getLastProp<0>()[1] = 0.0;
        }
        ++it;
    }
    Particles.map();
    Particles.ghost_get<0>();

    //We write the particles to check if the initialization is correct.
    Particles.write("Init");
    //! @cond [init1Subset] @endcond
    /**
     * @page Odeint_single_step Step by step time integration with Odeint
     *
     * ## Create the subset and Cast Global Pointers ## {#odeint_c_point}
     *
     *
     * On the particles we just created we need to constructed the subset object based on the numbering.
     * Further, We cast the Global Pointers so that Odeint RHS functor can recognize our openfpm distributed structure.
     *
     *
     * @snippet example/Numerics/OdeInt/main.cpp Pointer1Init
     */
    //! @cond [Pointer1Init] @endcond
    // Now we initialize the grid with a filled circle. Outside the circle, the value of Phi_0 will be -1, inside +1.
    //Now we construct the subsets based on the subset number.
    dist_vector_subset_type Particles_bulk(Particles, 0);
    dist_vector_subset_type Particles_boundary(Particles, 1);

    //We cast the global pointers to Particles and Particles_bulk as expected by the RHS functor.
    PointerDistGlobal = (void *) &Particles;
    PointerDistSubset = (void *) &Particles_bulk;
    //! @cond [Pointer1Init] @endcond

    /**
     *
     * @page Odeint_single_step Step by step time integration with Odeint
     *
     * ## Creating DCPSE Operators and aliases for expressions## {#odeint_c_dcpse}
     *
     * Here we create two dcpse based operators and alias the particle properties.
     *
     * @snippet example/Numerics/OdeInt/main.cpp DCPSE1Alias
     *
     */
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @cond [DCPSE1Alias] @endcond
    //We create the DCPSE Based operators for discretization of the operators.
    auto verletList = Particles.template getVerlet<VL_NON_SYMMETRIC|VL_SKIP_REF_PART>(rCut);
    Derivative_xx<decltype(verletList)> Dxx(Particles, verletList, 2, rCut);
    Derivative_yy<decltype(verletList)> Dyy(Particles, verletList, 2, rCut);
    //We create aliases for referring to property and and positions.
    auto Pos = getV<PROP_POS>(Particles);
    auto C = getV<0>(Particles);
    auto V = getV<1>(Particles);
    auto dC = getV<2>(Particles);

    //We create aliases for referring to the subset properties.
    auto C_bulk = getV<0>(Particles_bulk);
    auto V_bulk = getV<1>(Particles_bulk);
    auto dC_bulk = getV<2>(Particles_bulk);

    //! @cond [DCPSE1Alias] @endcond

    /**
     * @page Odeint_single_step Step by step time integration with Odeint
     *
     * ## Creating Odeint Objects ## {#odeint_c_obj}
     * Now we create a odeint stepper object (RK4 in this case. Please refer to odeint on more such methods). Since we are in 2d, we are going to use "state_type_2d_ofp". Which is a structure or state_type compatible with odeint. We further pass all the parameters including "boost::numeric::odeint::vector_space_algebra_ofp",which tell odeint to use openfpm algebra.
     * The template parameters are: state_type_2d_ofp (state type of X), double (type of the value inside the state), state_type_2d_ofp (state type of DxDt), double (type of the time), boost::numeric::odeint::vector_space_algebra_ofp (our algebra).
     *
     * We further create the an instance of the RHS Functor defined before main. This instance is needed by odeint to compute the stages.
     *
     * Also, we create the state type compatible with odeint and initialize the concentration in it.
     *
     * @snippet example/Numerics/OdeInt/main.cpp Odeint1I
     *
     */
    //! @cond [Odeint1I] @endcond
    //Now we create a odeint stepper object (RK4). Since we are in 2d, we are going to use "state_type_2d_ofp". Which is a structure or state_type compatible with odeint. We further pass all the parameters including "boost::numeric::odeint::vector_space_algebra_ofp",which tell odeint to use openfpm algebra.
    // The template parameters are: state_type_2d_ofp (state type of X), double (type of the value inside the state), state_type_2d_ofp (state type of DxDt), double (type of the time), boost::numeric::odeint::vector_space_algebra_ofp (our algebra)
    boost::numeric::odeint::runge_kutta4<state_type_2d_ofp, double, state_type_2d_ofp, double, boost::numeric::odeint::vector_space_algebra_ofp> Odeint_rk4;

    //The method Odeint_rk4 from Odeint, requires system (a function which computes RHS of the PDE), an instance of the Compute RHS functor. We create the System with the correct types and parameteres for the operators as declared before.
    RHSFunctor<Derivative_xx<decltype(verletList)>, Derivative_yy<decltype(verletList)>> System(Dxx, Dyy);

    //Furhter, odeint needs data in a state type "state_type_2d_ofp", we create one and fill in the initial condition.
    state_type_2d_ofp X;
    //Since we created a 2d state_type we initialize the two fields in the object data using the method get.
    X.data.get<x>() = C[x];
    X.data.get<y>() = C[y];
    //! @cond [Odeint1I] @endcond

    /**
    * @page Odeint_single_step Step by step time integration with Odeint
    *
    * ## Time Integration Loop ## {#odeint_c_2}
    * We create a counter for counting the step and initiliaze the time loop.
    *
    * Inside the time loop, we first compute the Velocity at the current position and observe the state of the sysem by writing the solution.
    *
    * We then Call the Odeint_rk4 method created above to do a rk4 step with arguments as the System, the state_type, current time t and the stepsize dt.
    *
    * Odeint updates X in place. hence we copy the computations of the time step back to Concentration (Bulk only as we dont want to change the boundary conditions).
    *
    * We then advect the particles (an Euler step) and do a map and ghost_get as needed after moving particles.
    *
    * We finally update the subset bulk and the DCPSE operators. And also update ctr and t.
    *
    * After the time loop. we deallocate the DCPSE operators and finalize the library.
    *
    * @snippet example/Numerics/OdeInt/main.cpp OdeintT
    *
    */
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @cond [OdeintT] @endcond
    // Now we will create a time loop for controlling the step size ourselves but call odeint to do the 4 stages of RK4.
    int ctr = 0;
    double t = 0, tf = 1, dt = 1e-2;
    while (t < tf) {
        //computing the velocity at the current position and at time step t (only in the bulk, so boundary remains 0).
        V_bulk[x] = -Pos[y] * exp(-10.0 * (Pos[x] * Pos[x] + Pos[y] * Pos[y]));
        V_bulk[y] = Pos[x] * exp(-10.0 * (Pos[x] * Pos[x] + Pos[y] * Pos[y]));
        //Observing the state
        Particles.deleteGhost();
        Particles.write_frame("PDE_Sol", ctr);
        Particles.ghost_get<0>();
        //calling rk4 with the function System, the state_type, current time t and the stepsize dt. It computes one step of the RK4.
        Odeint_rk4.do_step(System, X, t, dt);
        //Copying back the step, updating only the bulk values.
        C_bulk[x] = X.data.get<0>();
        C_bulk[y] = X.data.get<1>();
        //We advect the particles with the velocity using an Euler step.
        Pos = Pos + dt * V;
        //We need to map and get the ghost, as we moved the particles.
        Particles.map();
        Particles.ghost_get<0>();
        //We update the subset and operators as the particles moved.
        Particles_bulk.update();
        Particles_boundary.update();
        Particles.updateVerlet(verletList,rCut);
        Dxx.update(Particles);
        Dyy.update(Particles);
        //Reinitialzing as distributed size can change.
        X.data.get<0>()=C[x];
        X.data.get<1>()=C[y];
        ctr++;
        t += dt;
    } //time loop end

    //Deallocating the operators
    Dxx.deallocate(Particles);
    Dyy.deallocate(Particles);
    openfpm_finalize(); // Finalize openFPM library
    return 0;
} //main end
//! @cond [OdeintT] @endcond
/**
 * @page Odeint_single_step Step by step time integration with Odeint
 *
 * ## Full code ## {#odeint_c_full}
 *
 * @include example/Numerics/OdeInt/main.cpp
 */
