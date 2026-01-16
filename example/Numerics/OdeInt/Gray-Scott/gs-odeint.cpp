
// Include Vector Expression,Vector Expressions for Subset,DCPSE,Odeint header files
#include "Operators/Vector/vector_dist_operators.hpp"
#include "DCPSE/DCPSE_op/DCPSE_op.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"
#include <string>
//! @cond [Ode2Include] @endcond

constexpr int x = 0;
constexpr int y = 1;

double dt=1.0,tf=5000.0;

void *PointerDistGlobal;

typedef aggregate<VectorS<2, double>,VectorS<2, double>> Property_type;
typedef vector_dist<3, double, Property_type> dist_vector_type;


template<typename laplacian_type, typename verletList_type>
struct RHSFunctor
{

    //Intializing the operators
    laplacian_type &Lap;
    verletList_type &verletList;

    // Physical contants
    double K = 0.053;
    double F = 0.014;

    double d1 = 2*1e-4;
    double d2 = 1*1e-4;

    //Constructor
    RHSFunctor(laplacian_type &Lap, verletList_type& verletList) : Lap(Lap), verletList(verletList)
    {}

    void operator()( const state_type_2d_ofp &X , state_type_2d_ofp &dxdt , const double t ) const
    {
        //Casting the pointers to OpenFPM vector distributions
        dist_vector_type &Particles= *(dist_vector_type *) PointerDistGlobal;

        //Aliasing the properties.
        auto C = getV<0>(Particles);
        //These expressions only update the bulk values of C.
        C[x]=X.data.get<0>();
        C[y]=X.data.get<1>();
        Particles.ghost_get<0>(SKIP_LABELLING);
        // Particles.updateVerlet(verletList,verletList.getRCut());

        // We do the RHS computations for the Laplacian and reaction term
        // (Updating bulk only).
        dxdt.data.get<0>() = d1*Lap(C[x]) - C[x] * C[y] * C[y] + F - F * C[x];
        dxdt.data.get<1>() = d2*Lap(C[y]) + C[x] * C[y] * C[y] - (F+K) * C[y];
        //We copy back to the dxdt state_type for Odeint
        //=dC[x];
        //=dC[y];
    }
};

struct ObserverFunctor {

    int ctr;
    double t_old;

    //Constructor
    ObserverFunctor() {
        //a counter for counting the np. of steps
        ctr = 0;
        //Starting with t=0, we compute the step size take by t-t_old. So for the first observed step size is what we provide. Which will be 0-(-dt)=dt.
        t_old = -dt;
    }

    void operator()(state_type_2d_ofp &X,const double t) {
        if (ctr % 300 == 0) {
            dist_vector_type &Particles= *(dist_vector_type *) PointerDistGlobal;
            auto C = getV<0>(Particles);
            C[x]=X.data.get<0>();
            C[y]=X.data.get<1>();
            auto &v_cl=create_vcluster();
            if(v_cl.rank()==0)
            {
                std::cout<<"Time: "<<t<<", "<<"dt: "<<t-t_old<<std::endl;
            }
            Particles.deleteGhost(); 
            Particles.write_frame("PDE_sol",ctr,t);
            Particles.ghost_get<0>();
        }
    t_old=t;
    ctr++;        
    }

};


template <typename stepper_type, typename laplacian_type, typename verletList_type>
void run_stepper_const(dist_vector_type &Particles, std::vector<double> &runtime_v, laplacian_type &Lap, verletList_type& verletList) {

    RHSFunctor<laplacian_type, verletList_type> System(Lap, verletList);
    ObserverFunctor ObserveAndUpdate;
    auto C = getV<0>(Particles);
    auto InitC = getV<1>(Particles);

    state_type_2d_ofp x0;
    x0.data.get<x>() = InitC[x];
    x0.data.get<y>() = InitC[y];

    timer timer_integrate;
    timer_integrate.start();
    boost::numeric::odeint::integrate_const(stepper_type(), System, x0, 0.0, tf, dt, ObserveAndUpdate);
//    boost::numeric::odeint::integrate_const(stepper_type(), derivative, x0, 0.0, tf, dt);
    timer_integrate.stop();
    double rt=timer_integrate.getwct();
    auto &v_cl=create_vcluster();  
    v_cl.sum(rt);
    v_cl.execute();

    runtime_v.push_back(rt);
    C[x]=x0.data.get<x>();
    C[y]=x0.data.get<y>();
    if(v_cl.rank()==0)std::cout << "Runtime: " << rt << std::endl;
}

template <typename stepper_type, typename laplacian_type, typename verletList_type>
void run_stepper_adaptive(dist_vector_type &Particles, std::vector<double> &runtime_v, laplacian_type &Lap, verletList_type& verletList) {

    RHSFunctor<laplacian_type, verletList_type> System(Lap, verletList);
    ObserverFunctor ObserveAndUpdate;
    auto C = getV<0>(Particles);
    auto InitC = getV<1>(Particles);

    state_type_2d_ofp x0;
    x0.data.get<x>() = InitC[x];
    x0.data.get<y>() = InitC[y];

    timer timer_integrate;
    timer_integrate.start();
    boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(1e-3,1e-3,stepper_type()), System , x0 , 0.0 , tf , dt, ObserveAndUpdate);
//    boost::numeric::odeint::integrate_adaptive( stepper_type() , derivative , x0 , 0.0 , tf , dt);
    timer_integrate.stop();
    double rt=timer_integrate.getwct();
    auto &v_cl=create_vcluster();  
    v_cl.sum(rt);
    v_cl.execute();
    runtime_v.push_back(rt);

    C[x]=x0.data.get<x>();
    C[y]=x0.data.get<y>();
    if(v_cl.rank()==0)std::cout << "Runtime: " << rt << std::endl;
}

double average(std::vector<double> &nums) {
    return std::accumulate(nums.begin(), nums.end(), 0.0) / static_cast<double>(nums.size());
}


int main(int argc, char *argv[])
{
    //	initialize library
    openfpm_init(&argc, &argv);
    tf=std::atof(argv[2]);
    // output
    std::vector<double> runtime_rk4_const;
    std::vector<double> runtime_rk5_const;
    std::vector<double> runtime_rk78_const;
    std::vector<double> runtime_rk5_adapt;
    size_t gdsz=std::atof(argv[1]);
    Box<3,double> box({0.0,0.0,0.0},{2.5,2.5,2.5});
    size_t sz[3] = {gdsz,gdsz,gdsz};
    // Define periodicity of the grid
    size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};
    double spacing[3];
    spacing[0] = 2.5 / (sz[0]);
    spacing[1] = 2.5 / (sz[1]);
    spacing[2] = 2.5 / (sz[2]);
    double rCut = 2.9 * spacing[0];
    Ghost<3, double> ghost(rCut);

    dist_vector_type Particles(0, box, bc, ghost);
    Particles.setPropNames({"Concentration","Initial"});

    auto it = Particles.getGridIterator(sz);
    while (it.isNext()) {
        Particles.add();
        auto key = it.get();
        double x = 0.0 + key.get(0) * spacing[0];
        Particles.getLastPos()[0] = x;
        double y = 0.0 + key.get(1) * spacing[1];
        Particles.getLastPos()[1] = y;
        double z = 0.0 + key.get(2) * spacing[2];
        Particles.getLastPos()[2] = z;
        // Here fill the Initial value of the concentration.
        Particles.template getLastProp<1>()[0] = 1.0;
        Particles.template getLastProp<1>()[1] = 0.0;

        if (x > 1.55 && x < 1.85 && y > 1.55 && y < 1.85 && z > 1.55 && z < 1.85) {
            Particles.template getLastProp<1>()[0] = 0.5 + (((double)std::rand())/RAND_MAX -0.5)/10.0;
            Particles.template getLastProp<1>()[1] = 0.25 + (((double)std::rand())/RAND_MAX -0.5)/20.0;
        }

        ++it;
    }
    Particles.map();
    Particles.ghost_get<0>();


    // Now we initialize the grid with a filled circle. Outside the circle, the value of Phi_0 will be -1, inside +1.
    //Now we construct the subsets based on the subset number.

    //We cast the global pointers to Particles and Particles_bulk as expected by the RHS functor.
    PointerDistGlobal = (void *) &Particles;
    auto verletList = Particles.template getVerlet<VL_NON_SYMMETRIC|VL_SKIP_REF_PART>(rCut);
    //We create the DCPSE Based Laplacian operator.

    Laplacian<decltype(verletList)> Lap(Particles, verletList, 2, rCut, support_options::RADIUS);
    auto C = getV<0>(Particles);
    auto Init = getV<1>(Particles);
    C=Init;
    //Now we create a odeint stepper object (RK4). Since we are in 2d, we are going to use "state_type_2d_ofp". Which is a structure or state_type compatible with odeint. We further pass all the parameters including "boost::numeric::odeint::vector_space_algebra_ofp",which tell odeint to use openfpm algebra.
    // The template parameters are: state_type_2d_ofp (state type of X), double (type of the value inside the state), state_type_2d_ofp (state type of DxDt), double (type of the time), boost::numeric::odeint::vector_space_algebra_ofp (our algebra)
    typedef boost::numeric::odeint::runge_kutta4<state_type_2d_ofp, double, state_type_2d_ofp, double, boost::numeric::odeint::vector_space_algebra_ofp> Odeint_rk4;
    typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_type_2d_ofp,double,state_type_2d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp> Odeint_rk5;
    typedef boost::numeric::odeint::runge_kutta_fehlberg78< state_type_2d_ofp,double,state_type_2d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp> Odeint_rk8;
    //The method Odeint_rk4 from Odeint, requires system (a function which computes RHS of the PDE), an instance of the Compute RHS functor. We create the System with the correct types and parameteres for the operators as declared before.
    RHSFunctor<Laplacian<decltype(verletList)>, decltype(verletList)> System(Lap, verletList);

    //Since we are using Odeint to control the time steps, we also create a observer instance. Which also updates the position via an euler step for moving thr particles.
    ObserverFunctor ObserveAndUpdate;


    //Furhter, odeint needs data in a state type "state_type_2d_ofp", we create one and fill in the initial condition.
    state_type_2d_ofp X;
    //Since we created a 2d state_type we initialize the two fields in the object data using the method get.
    X.data.get<x>() = C[0];
    X.data.get<y>() = C[1];


    std::vector<double> inter_times; // vector to store intermediate time steps taken by odeint.
        Particles.deleteGhost();
        Particles.write("Initial");
        Particles.ghost_get<0>();
    //for (int i = 0; i < 3; ++i) {
        // run_stepper_const<Odeint_rk4>(Particles, runtime_rk4_const,Lap, verletList);
        // Particles.deleteGhost();
        // Particles.write("RK4final");
        // Particles.ghost_get<0>();
        // run_stepper_const<Odeint_rk5>(Particles, runtime_rk5_const,Lap,verletList);
        // Particles.deleteGhost();
        // Particles.write("RK5final");
        // Particles.ghost_get<0>();
        // run_stepper_const<Odeint_rk8>(Particles, runtime_rk78_const,Lap,verletList);
        // Particles.deleteGhost();
        // Particles.write("RK8");
        // Particles.ghost_get<0>();
        run_stepper_adaptive<Odeint_rk5>(Particles, runtime_rk5_adapt,Lap,verletList);
        Particles.deleteGhost();
        Particles.write("AdapRK5");
        Particles.ghost_get<0>();
        


    //}

//    size_t steps = boost::numeric::odeint::integrate_const(Odeint_rk4, System, X, 0.0, tf, dt, ObserveAndUpdate);
//    size_t steps = boost::numeric::odeint::integrate_adaptive( boost::numeric::odeint::make_controlled( 1.0e-7 , 1.0e-7 , Odeint_rk5()) , System , X , 0.0 , tf , dt, ObserveAndUpdate );
    
    auto & vcl = create_vcluster();

    if (vcl.getProcessUnitID() == 0) {

        std::ofstream output_runtime;

        // head line
        std::string headline = "cores,rk4,dopri5,fehlberg78,dopri5 adaptive";
        if (access("runtime.csv", F_OK) == -1) {
            output_runtime.open("runtime.csv");
            output_runtime << headline << "\n";
        } else {
            output_runtime.open("runtime.csv", std::ios::app);
        }

        output_runtime << vcl.getProcessingUnits()
                       << "," << average(runtime_rk4_const)
                       << "," << average(runtime_rk5_const)
                       << "," << average(runtime_rk78_const)
                       << "," << average(runtime_rk5_adapt)
                       << "\n";
    }

    //Deallocating the operators
    Lap.deallocate(Particles);
    openfpm_finalize(); // Finalize openFPM library
    return 0;
} //main end
