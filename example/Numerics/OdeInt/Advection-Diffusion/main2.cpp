//
// Created by Abhinav Singh @absingh on 26/04/2021
//
/**
 * @page Odeint_multiple_steps Multiple steps integration with Odeint
 * *
 * In this example, we perform time integration in a 2d domain of particles of a following partial differential equation:
 *
 * @f[ \frac{\partial\vec{C}}{dt}=\vec{V}.\nabla_{\{x,y\}} U + 0.1*\Delta_{\{x,y\}} U @f]
 * in 2d domain [-1,-1]*[1,1] with a fixed velocity \f$\vec{V}(x,y)=(-ye^{-10(x^2+y^2)},xe^{-10(x^2+y^2)})\f$, and the boundary conditions on the walls as no-slip for velocity \f$\vec{V}=0\f$ and sink for the chemicals \f$\vec{C}=0\f$ for all time \f$t\f$.
 *  Further, we start with the initial condition for the concentration as
 *   @f[\vec{C}=\begin{cases}
 *   (1,0)\text{ for } x=0,-0.5<y<0\\
 *   (0,1)\text{ for } x=0, 0<y<0.5\\
 *   (0,0)  \text{ for the rest of the domain}\\
 * \end{cases} @f]
 *
 * We do that by emplying a Lagrangian frame of reference. Hence the Pde is transformed to:
 *
 * @f[\begin{align}
 *   \frac{\partial\vec{X}}{dt}=\vec{V}\\
 *   \frac{\partial\vec{C}}{dt}=0.1*\Delta_{\{x,y\}} U
 * \end{align} @f]
 * This is a system of PDEs. We decouple the moving of the particles from the evolution of the chemicals as it can be expensive to recompute derivatives at every stage of a time integrator in a single step.
 *
 * Output:
 * Time series data of the PDE Solution.
 *
 */

/**
 * @page Odeint_multiple_steps Multiple steps integration with Odeint
 *
 * ## Includes ## {#ode_c2_include}
 *
 * These are the header files that we need to include:
 *
 * @snippet example/Numerics/Odeint/main2.cpp Ode2Include
 *
 */
//! @cond [Ode2Include] @endcond

// Include Vector Expression,Vector Expressions for Subset,DCPSE,Odeint header files
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "DCPSE/DCPSE_op/DCPSE_op.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"
//! @cond [Ode2Include] @endcond


/**
 * @page Odeint_multiple_steps Multiple steps integration with Odeint
 *
 * ## Initialization of the global parameters## {#ode_c2_init}
 *
 * We start with
 * * Initializing certain global parameteres we will use: such as x,y to refer to the dimensions 0 and 1. (Makes it easier to read equations in code)
 * * Creating empty pointers for coupling openfpm distributed vector with odeint. One for the entire distributed vector and another for the subset or bulk.
 * We seperate bulk and the entire distribution as it makes it easier to impose boundary conditions. (Which will be more apparant in ComputeRHS of the PDE)
 * Note that a subset expression always comes at the left hand side of a computation. (The semantics of the expressions is by denoting what we want to update from regular expressions)
 * Creating aliases of the types of the datasructures we are going to use in OpenFPM. (Bulk Alias can only be used on LHS, say to update onl the bulk)
 *
 * Property_type as the type of properties we wish to use.
 *
 * dist_vector_type as the 2d openfpm distributed vector type
 *
 * dist_vector_type as the 2d openfpm distributed subset vector type
 *
 * @snippet example/Numerics/Odeint/main2.cpp Initialization__two
 *
 */
//! @cond [Initialization__two] @endcond

constexpr int x = 0;
constexpr int y = 1;

double dt=1e-2,tf=1.0,vf=1.0;

void *PointerDistGlobal, *PointerDistSubset,*PointerDistSubset2;

typedef aggregate<VectorS<2, double>, VectorS<2, double>, VectorS<2, double>> Property_type;
typedef vector_dist_ws<2, double, Property_type> dist_vector_type;
typedef vector_dist_subset<2, double, Property_type> dist_vector_subset_type;
//! @cond [Initialization__two] @endcond

/**
 * @page Odeint_multiple_steps Multiple steps integration with Odeint
 *
 * ## Creating the RHS Functor## {#ode_c2_rhs}
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
 * We pas on the openfpm distributed state types as
 * void operator()( const state_type_2d_ofp &X , state_type_2d_ofp &dxdt , const double t ) const
 *
 * Since we would like to use openfpm here. We cast back the global pointers created before to access the Openfpm distributed vector here.
 * (Note that these pointers needs to initialized in the main(). Further, 'state_type_2d_ofp' is a temporal structure, which means it does not have the ghost.
 * Hence we copy the current state back to the openfpm vector from the openfpm state type X.
 * We do our computations as required.
 * Then we copy back the output into the state_type dxdt
 *
 * )
 *
 * @snippet example/Numerics/Odeint/main2.cpp RHS2Functor
 *
 */
//! @cond [RHS2Functor] @endcond
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
//! @cond [RHS2Functor] @endcond


/**
 * @page Odeint_multiple_steps Multiple steps integration with Odeint
 *
 * ## Creating the Observer Functor## {#ode_c2_obs}
 *
 * There are multiple ways in which the system can be integrated. For example, and ideally, we could put both the PDEs into the RHS functor (Moving the particles at every stage). This can be expensive.
 * However, Often in numerical simulations, Both the PDEs can be integrated with seperate steppers.
 * To achieve this we will use the observer functor.
 * The observer functor is called before every time step evolution by odeint. Hence we can use it to update the position of the particles, with an euler step. and also update the operators and write/observe the solution.
 *
 * Now we create the Observer functor. Please refer to ODE_int for more detials.
 * Note that we have templated it with two operator types DXX and DYY again. But we never use them. It is just an example to show how versatile the observer can be.
 *
 * All Observer computations  needs to happen in the operator ().
 * Odeint expects the arguments here to be an input state_type X, and time t.
 * We pass on the openfpm distributed state types as
 * void operator()( const state_type_2d_ofp &X , const double t ) const
 *
 * Since we would like to use again openfpm here. We cast back the global pointers created before to access the Openfpm distributed vector here.
 * (Note that these pointers needs to initialized in the main(). Further, 'state_type_2d_ofp' is a temporal structure, which means it does not have the ghost.
 * Hence we copy the current state back to the openfpm vector from the openfpm state type X.
 * We do our computations as required.
 * Then we copy back the output into the state_type dxdt.
 *
 * @snippet example/Numerics/Odeint/main2.cpp Observer2Functor
 *
 */
//! @cond [Observer2Functor] @endcond
template<typename DXX, typename DYY, typename VerletList_type>
struct ObserverFunctor {

    DXX &Dxx;
    DYY &Dyy;
    VerletList_type &verletList;
    int ctr;
    double t_old;
    double rCut;

    //Constructor
    ObserverFunctor(DXX &Dxx, DYY &Dyy, VerletList_type& verletList, double rCut) : Dxx(Dxx), Dyy(Dyy), verletList(verletList), rCut(rCut) {
        //a counter for counting the np. of steps
        ctr = 0;
        //Starting with t=0, we compute the step size take by t-t_old. So for the first observed step size is what we provide. Which will be 0-(-dt)=dt.
        t_old = -dt;
    }

    void operator()(state_type_2d_ofp &X, double t) {
        //Casting the pointers to OpenFPM vector distributions
        dist_vector_type &Particles = *(dist_vector_type *) PointerDistGlobal;
        dist_vector_subset_type &Particles_bulk = *(dist_vector_subset_type *) PointerDistSubset;
        dist_vector_subset_type &Particles_boundary = *(dist_vector_subset_type *) PointerDistSubset2;

        //Aliasing the position and properties.
        auto Pos = getV<PROP_POS>(Particles);
        auto Concentration = getV<0>(Particles);
        auto Velocity = getV<1>(Particles);
        auto Concentration_bulk = getV<0>(Particles_bulk);
        auto Velocity_bulk = getV<1>(Particles_bulk);

        //We would like to move the particles after t=0. (Remember, odeint calls the observer before taking the step.)
        if (t != 0) {
            //copy back the state after the time step into data structure. This is required as we are going to move the particles and the distributed state can be resized correctly (by copying back after map). Also this expression preserves the boundary condition on concentration.
            Concentration_bulk[x] = X.data.get<0>();
            Concentration_bulk[y] = X.data.get<1>();
            //computing the velocity and move the particles
            Velocity_bulk[x] = -vf*Pos[y] * exp(-10.0 * (Pos[x] * Pos[x] + Pos[y] * Pos[y]));
            Velocity_bulk[y] = vf*Pos[x] * exp(-10.0 * (Pos[x] * Pos[x] + Pos[y] * Pos[y]));
            Pos = Pos + dt * Velocity;
            //Map and ghost_get is required after moving particles.
            Particles.map();
            Particles.ghost_get<0>();
            //Updating the subset and operators based on new positions
            Particles_bulk.update();
            Particles_boundary.update();
            Particles.updateVerlet(verletList,rCut);
            Dxx.update(Particles);
            Dyy.update(Particles);

            //Since we did Map, we assign the Odeint state again.
            X.data.get<0>() = Concentration[x];
            X.data.get<1>() = Concentration[y];
            //counting the step number

        }
        ctr++;
        std::cout<<"Taking a step at t="<<t<<" with dt="<<t-t_old<<std::endl;
        t_old=t;
        //Writing the data after every step
        Particles.deleteGhost();
        Particles.write_frame("PDE_sol", ctr);
        Particles.ghost_get<0>();

    }
};
//! @cond [Observer2Functor] @endcond
/**
 * @page Odeint_multiple_steps Multiple steps integration with Odeint
 *
 * ## Initializating OpenFPM ## {#odeint_c2_initmain}
 *
 * We start with
 * * Initializing OpenFPM
 *
 * @snippet example/Numerics/Odeint/main2.cpp initParticles2
 *
 */
//! @cond [initParticles2] @endcond
int main(int argc, char *argv[])
{
    //	initialize library
    openfpm_init(&argc, &argv);

    dt=std::atof(argv[1]);
    tf=std::atof(argv[2]);
    vf=std::atof(argv[3]);
    //! @cond [initParticles2] @endcond

    /**
     * @page Odeint_multiple_steps Multiple steps integration with Odeint
     *
     * ## Creating Particles and assigning subsets ## {#odeint_c2_indices}
     *
     * We create a particle distribution we certain rCut for the domain [-1,-1] to [1,1].
     *
     * Also, we fill the initial concentration as C_1(x=0,y>0 & y<0.5,t=0)=1,C_2(x=0,y<0 & y>-0.5,t=0)=1 and 0 everywhere else.
     *
     * @snippet example/Numerics/Odeint/main2.cpp init2Subset
     *
     */
    //! @cond [init2Subset] @endcond
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
    //! @cond [init2Subset] @endcond
    /**
     * @page Odeint_multiple_steps Multiple steps integration with Odeint
     *
     * ## Create the subset and Cast Global Pointers ## {#odeint_c2_point}
     *
     *
     * On the particles we just created we need to constructed the subset object based on the numbering.
     * Further, We cast the Global Pointers so that Odeint RHS functor can recognize our openfpm distributed structure.
     *
     *
     * @snippet example/Numerics/Odeint/main2.cpp Pointer2Init
     */
    //! @cond [Pointer2Init] @endcond
    // Now we initialize the grid with a filled circle. Outside the circle, the value of Phi_0 will be -1, inside +1.
    //Now we construct the subsets based on the subset number.
    dist_vector_subset_type Particles_bulk(Particles, 0);
    dist_vector_subset_type Particles_boundary(Particles, 1);

    //We cast the global pointers to Particles and Particles_bulk as expected by the RHS functor.
    PointerDistGlobal = (void *) &Particles;
    PointerDistSubset = (void *) &Particles_bulk;
    PointerDistSubset2 = (void *) &Particles_boundary;
    //! @cond [Pointer2Init] @endcond

    /**
     * @page Odeint_multiple_steps Multiple steps integration with Odeintx`
     *
     * ## Creating DCPSE Operators and aliases for expressions## {#odeint_c2_dcpse}
     *
     * Here we create two dcpse based operators and alias the particle properties.
     *
     * @snippet example/Numerics/Odeint/main2.cpp DCPSE2Alias
     *
     */
    //! @cond [DCPSE2Alias] @endcond
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
    auto dC_bulk = getV<2>(Particles_bulk);
    //! @cond [DCPSE2Alias] @endcond

    /**
     * @page Odeint_multiple_steps Multiple steps integration with Odeintx`
     *
     * ## Creating Odeint Objects## {#odeint_c2_obj}
     * Now we create a odeint stepper object (RK4 in this case. Please refer to odeint on more such methods or examples listed as comments after calling the method). Since we are in 2d, we are going to use "state_type_2d_ofp". Which is a structure or state_type compatible with odeint. We further pass all the parameters including "boost::numeric::odeint::vector_space_algebra_ofp",which tell odeint to use openfpm algebra.
     * The template parameters are: state_type_2d_ofp (state type of X), double (type of the value inside the state), state_type_2d_ofp (state type of DxDt), double (type of the time), boost::numeric::odeint::vector_space_algebra_ofp (our algebra).
     *
     * We further create the an instance of the RHS Functor defined before main. This instance is needed by odeint to compute the stages.
     *
     * Also, we create the state type compatible with odeint and initialize the concentration in it.
     *
     * @snippet example/Numerics/Odeint/main2.cpp Odeint2I
     *
     */
    //! @cond [Odeint2I] @endcond
    //Now we create a odeint stepper object (RK4). Since we are in 2d, we are going to use "state_type_2d_ofp". Which is a structure or state_type compatible with odeint. We further pass all the parameters including "boost::numeric::odeint::vector_space_algebra_ofp",which tell odeint to use openfpm algebra.
    // The template parameters are: state_type_2d_ofp (state type of X), double (type of the value inside the state), state_type_2d_ofp (state type of DxDt), double (type of the time), boost::numeric::odeint::vector_space_algebra_ofp (our algebra)
    boost::numeric::odeint::runge_kutta4<state_type_2d_ofp, double, state_type_2d_ofp, double, boost::numeric::odeint::vector_space_algebra_ofp> Odeint_rk4;

    //The method Odeint_rk4 from Odeint, requires system (a function which computes RHS of the PDE), an instance of the Compute RHS functor. We create the System with the correct types and parameteres for the operators as declared before.
    RHSFunctor<Derivative_xx<decltype(verletList)>, Derivative_yy<decltype(verletList)>> System(Dxx, Dyy);

    //Since we are using Odeint to control the time steps, we also create a observer instance. Which also updates the position via an euler step for moving thr particles.
    ObserverFunctor<Derivative_xx<decltype(verletList)>, Derivative_yy<decltype(verletList)>, decltype(verletList)> ObserveAndUpdate(Dxx, Dyy, verletList, rCut);


    //Furhter, odeint needs data in a state type "state_type_2d_ofp", we create one and fill in the initial condition.
    state_type_2d_ofp X;
    //Since we created a 2d state_type we initialize the two fields in the object data using the method get.
    X.data.get<x>() = C[0];
    X.data.get<y>() = C[1];
    //! @cond [Odeint2I] @endcond

    /**
    * @page Odeint_multiple_steps Multiple steps integration with Odeintx`
    *
    * ## Calling Odeint ## {#odeint_c2_1211}
    * We initiliaze the time variable t, step_size dt and final time tf.
    *
    * We create a vector for storing the intermidiate time steps, as most odeint calls return such an object.
    *
    * We then Call the Odeint_rk4 method created above to do a rk4 time integration from t0 to tf with arguments as the System, the state_type, current time t and the stepsize dt.
    *
    * Odeint updates X in place. And automatically advect the particles (an Euler step) and do a map and ghost_get as needed after moving particles by calling the observer.
    *
    * The observer also update the subset bulk and the DCPSE operators.
    *
    * We finally deallocate the DCPSE operators and finalize the library.
    *
    * @snippet example/Numerics/Odeint/main2.cpp OdeintTCall
    *
    */
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @cond [OdeintTCall] @endcond
    std::vector<double> inter_times; // vector to store intermediate time steps taken by odeint.

    size_t steps = boost::numeric::odeint::integrate_const(Odeint_rk4, System, X, 0.0, tf, dt, ObserveAndUpdate);
    /* From Odeint: (Please refer to Odeint on more details about these steppers.)
    runge_kutta_dopri5 is maybe the best default stepper. It has step size control as well as dense-output functionality. Simple create a dense-output stepper by make_dense_output( 1.0e-6 , 1.0e-5 , runge_kutta_dopri5< state_type >() ).
    runge_kutta4 is a good stepper for constant step sizes. It is widely used and very well known. If you need to create artificial time series this stepper should be the first choice.
    'runge_kutta_fehlberg78' is similar to the 'runge_kutta4' with the advantage that it has higher precision. It can also be used with step size control.
    adams_bashforth_moulton is very well suited for ODEs where the r.h.s. is expensive (in terms of computation time). It will calculate the system function only once during each step.*/
    /*Below are listed examples of using different steppers and their usage.*/
    //size_t steps = integrate_adaptive(controlled_stepper_type() , System , X , t , tf , dt, ObserveAndUpdate);
    //size_t steps = boost::numeric::odeint::integrate_adaptive( boost::numeric::odeint::make_controlled( 1.0e-7 , 1.0e-7 ,boost::numeric::odeint::runge_kutta_cash_karp54< state_type_2d_ofp,double,state_type_2d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp>() ) , System , X , t , tf , dt, ObserveAndUpdate );
    //size_t steps = boost::numeric::odeint::integrate_const( boost::numeric::odeint::make_dense_output(1.0e-7, 1.0e-7,  boost::numeric::odeint::runge_kutta_dopri5< state_type_2d_ofp,double,state_type_2d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp >() ) , System , X , t , tf , dt, ObserveAndUpdate );
    //size_t steps = boost::numeric::odeint::integrate_adaptive( boost::numeric::odeint::make_dense_output( 1.0e-8 , 1.0e-8 ,  boost::numeric::odeint::runge_kutta_dopri5< state_type_2d_ofp,double,state_type_2d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp >() ) , System , X , t , tf , dt, ObserveAndUpdate );
    //size_t steps = boost::numeric::odeint::integrate_adaptive( boost::numeric::odeint::make_controlled( 1.0e-7 , 1.0e-7 ,  boost::numeric::odeint::runge_kutta_dopri5< state_type_2d_ofp,double,state_type_2d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp >() ) , System , X , t , tf , dt, ObserveAndUpdate );


    std::cout << "No. of Time steps taken: " << steps << std::endl;
    //Copying back the final solution and outputting the number of steps taken by Odeint.
    C_bulk[x] = X.data.get<0>();
    C_bulk[y] = X.data.get<1>();
    //Deallocating the operators
    Dxx.deallocate(Particles);
    Dyy.deallocate(Particles);
    openfpm_finalize(); // Finalize openFPM library
    return 0;
} //main end
//! @cond [OdeintTCall] @endcond
/**
 * @page Odeint_multiple_steps Multiple steps integration with Odeint
 *
 * ## Full code ## {#odeint_c2_full}
 *
 * @include example/Numerics/Odeint/main2.cpp
 */