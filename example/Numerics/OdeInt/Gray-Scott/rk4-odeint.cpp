
#include <stddef.h>
#include "Vector/vector_dist.hpp"
#include <iostream>
#include <fstream>
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "Decomposition/Distribution/SpaceDistribution.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"
#include <vector>
#include <numeric>

double t=-5.0, tf=5.0, dt = 0.01;


void Exponential2( const state_type_1d_ofp &x , state_type_1d_ofp &dxdt , const double t )
{
    dxdt.data.get<0>() = x.data.get<0>();
}

void sigmoid2( const state_type_1d_ofp &x , state_type_1d_ofp &dxdt , const double t )
{
    // g'(t) = g(t) * (1-sigma)
    dxdt.data.get<0>() = x.data.get<0>() * (1- (1/(1+exp(-t))) );
}


template <typename stepper_type>
void run_stepper_const(vector_dist<2, double, aggregate<double, double,double>> &Particles,
                 void (*derivative)(const state_type_1d_ofp &, state_type_1d_ofp &, const double),
                 std::vector<double> &runtime_v,std::vector<double> &norm_inf_v, std::vector<double> &norm_l2_v) {

    auto Init = getV<0>(Particles);
    auto Sol = getV<1>(Particles);
    auto OdeSol = getV<2>(Particles);

    state_type_1d_ofp x0;
    x0.data.get<0>() = Init;

    timer timer_integrate;
    timer_integrate.start();
    boost::numeric::odeint::integrate_const(stepper_type(), derivative, x0, t, tf, dt);
    timer_integrate.stop();


    OdeSol = x0.data.get<0>();

    double norm_inf = 0.0;
    double norm_l2 = 0.0;

    auto it2 = Particles.getDomainIterator();
    while (it2.isNext()) {
        auto p = it2.get();
        // calculate error
        double diff = Particles.getProp<1>(p) - Particles.getProp<2>(p);
        // calculate infinity norm
        if (fabs(diff) > norm_inf) {
            norm_inf = fabs(diff);
        }
        // calculate L2 norm
        norm_l2 += diff * diff;
        ++it2;
    }
    norm_l2 = sqrt(norm_l2);
    double rt=timer_integrate.getwct();
    auto &v_cl=create_vcluster();  
    v_cl.sum(rt);
    v_cl.execute();

    runtime_v.push_back(rt);
    norm_inf_v.push_back(norm_inf);
    norm_l2_v.push_back(norm_l2);

    //std::cout << "Runtime: " << rt/v_cl.size() << std::endl;
    //std::cout << "Inf: " << norm_inf << std::endl;
    //std::cout << "L2: " << norm_l2 << std::endl;

}

template <typename stepper_type>
void run_stepper_adaptive(vector_dist<2, double, aggregate<double, double,double>> &Particles,
                       void (*derivative)(const state_type_1d_ofp &, state_type_1d_ofp &, const double),
                       std::vector<double> &runtime_v,std::vector<double> &norm_inf_v, std::vector<double> &norm_l2_v) {

    auto Init = getV<0>(Particles);
    auto Sol = getV<1>(Particles);
    auto OdeSol = getV<2>(Particles);

    state_type_1d_ofp x0;
    x0.data.get<0>() = Init;

    timer timer_integrate;
    timer_integrate.start();
    boost::numeric::odeint::integrate_adaptive( stepper_type() , derivative , x0 , t , tf , dt);
    timer_integrate.stop();


    OdeSol = x0.data.get<0>();

    double norm_inf = 0.0;
    double norm_l2 = 0.0;

    auto it2 = Particles.getDomainIterator();
    while (it2.isNext()) {
        auto p = it2.get();
        // calculate error
        double diff = Particles.getProp<1>(p) - Particles.getProp<2>(p);
        // calculate infinity norm
        if (fabs(diff) > norm_inf) {
            norm_inf = fabs(diff);
        }
        // calculate L2 norm
        norm_l2 += diff * diff;
        ++it2;
    }
    norm_l2 = sqrt(norm_l2);

    double rt=timer_integrate.getwct();
    auto &v_cl=create_vcluster();  
    v_cl.sum(rt);
    v_cl.execute();

    runtime_v.push_back(rt);
    norm_inf_v.push_back(norm_inf);
    norm_l2_v.push_back(norm_l2);

    //std::cout << "Runtime: " << rt/v_cl.size() << std::endl;
    //std::cout << "Inf: " << norm_inf << std::endl;
    //std::cout << "L2: " << norm_l2 << std::endl;

}

double average(std::vector<double> &nums) {
    return std::accumulate(nums.begin(), nums.end(), 0.0) / static_cast<double>(nums.size());
}

int main(int argc, char* argv[])
{

    // initialize the library
    openfpm_init(&argc,&argv);

    // std::cout << "odeint"<< std::endl;


    // output
    std::vector<double> runtime_rk4_const;
    std::vector<double> runtime_rk5_const;
    std::vector<double> runtime_rk78_const;
    std::vector<double> runtime_rk5_adapt;

    std::vector<double> norm_inf_rk4_const;
    std::vector<double> norm_inf_rk5_const;
    std::vector<double> norm_inf_rk78_const;
    std::vector<double> norm_inf_rk5_adapt;

    std::vector<double> norm_l2_rk4_const;
    std::vector<double> norm_l2_rk5_const;
    std::vector<double> norm_l2_rk78_const;
    std::vector<double> norm_l2_rk5_adapt;

    size_t edgeSemiSize = std::stof(argv[1]);
    const size_t sz[2] = {edgeSemiSize,edgeSemiSize };
    Box<2, double> box({ 0, 0 }, { 1.0, 1.0 });
    size_t bc[2] = { NON_PERIODIC, NON_PERIODIC };
    double spacing[2];
    spacing[0] = 1.0 / (sz[0] - 1);
    spacing[1] = 1.0 / (sz[1] - 1);
    double rCut = 3.9 * spacing[0];
    Ghost<2, double> ghost(rCut);

    vector_dist<2, double, aggregate<double, double,double>> Particles(0, box, bc, ghost);

    // analytical solution
    auto it = Particles.getGridIterator(sz);
    while (it.isNext()) {
        Particles.add();
        auto key = it.get();
        mem_id k0 = key.get(0);
        double xp0 = k0 * spacing[0];
        Particles.getLastPos()[0] = xp0;
        mem_id k1 = key.get(1);
        double yp0 = k1 * spacing[1];
        Particles.getLastPos()[1] = yp0;

        // Particles.getLastProp<0>() = xp0*yp0*1/(1+exp(-t));
        // Particles.getLastProp<1>() = xp0*yp0*1/(1+exp(-tf));

        Particles.getLastProp<0>() = xp0 * yp0 * exp(t);
        Particles.getLastProp<1>() = xp0 * yp0 * exp(tf);

        ++it;
    }

    auto Init = getV<0>(Particles);
    auto Sol = getV<1>(Particles);
    auto OdeSol = getV<2>(Particles);

    state_type_1d_ofp x0;

    typedef boost::numeric::odeint::runge_kutta4<state_type_1d_ofp, double, state_type_1d_ofp, double, boost::numeric::odeint::vector_space_algebra_ofp> rk4_stepper_type;
    typedef boost::numeric::odeint::runge_kutta_dopri5<state_type_1d_ofp, double, state_type_1d_ofp, double, boost::numeric::odeint::vector_space_algebra_ofp> dopri5_stepper_type;
    typedef boost::numeric::odeint::runge_kutta_fehlberg78<state_type_1d_ofp, double, state_type_1d_ofp, double, boost::numeric::odeint::vector_space_algebra_ofp> fehlberg78_stepper_type;

//    void (*derivative)(const state_type_1d_ofp &, state_type_1d_ofp &, const double);
//    derivative = &Exponential2;

    for (int i = 0; i < 3; ++i) {

        run_stepper_const<rk4_stepper_type>(Particles, &Exponential2, runtime_rk4_const, norm_inf_rk4_const, norm_l2_rk4_const);
        run_stepper_const<dopri5_stepper_type>(Particles, &Exponential2, runtime_rk5_const, norm_inf_rk5_const, norm_l2_rk5_const);
        run_stepper_const<fehlberg78_stepper_type>(Particles, &Exponential2, runtime_rk78_const, norm_inf_rk78_const, norm_l2_rk78_const);
        run_stepper_adaptive<dopri5_stepper_type>(Particles, &Exponential2, runtime_rk5_adapt, norm_inf_rk5_adapt, norm_l2_rk5_adapt);

    }

    auto & vcl = create_vcluster();

    if (vcl.getProcessUnitID() == 0) {


        std::ofstream output_runtime;
        std::ofstream output_inf_norm;
        std::ofstream output_l2_norm;


        // head line
        std::string headline = "cores,rk4,dopri5,fehlberg78,dopri5 adaptive";
        if (access( "runtime.csv", F_OK ) == -1) {
            output_runtime.open("runtime.csv");
            output_runtime << headline << "\n";
        }
        else {
            output_runtime.open("runtime.csv",std::ios::app);
        }

        if (access( "inf_norm.csv", F_OK ) == -1) {
            output_inf_norm.open("inf_norm.csv");
            output_inf_norm << headline << "\n";
        }
        else {
            output_inf_norm.open("inf_norm.csv",std::ios::app);
        }

        if (access( "l2_norm.csv", F_OK ) == -1) {
            output_l2_norm.open("l2_norm.csv");
            output_l2_norm << headline << "\n";
        }
        else {
            output_l2_norm.open("l2_norm.csv",std::ios::app);
        }

        output_runtime << vcl.getProcessingUnits()
            << "," << average(runtime_rk4_const)
            << "," << average(runtime_rk5_const)
            << "," << average(runtime_rk78_const)
            << "," << average(runtime_rk5_adapt)
            << "\n";

        output_inf_norm << vcl.getProcessingUnits()
            << "," << average(norm_inf_rk4_const)
            << "," << average(norm_inf_rk5_const)
            << "," << average(norm_inf_rk78_const)
            << "," << average(norm_inf_rk5_adapt)
            << "\n";

        output_l2_norm << vcl.getProcessingUnits()
            << "," << average(norm_l2_rk4_const)
            << "," << average(norm_l2_rk5_const)
            << "," << average(norm_l2_rk78_const)
            << "," << average(norm_l2_rk5_adapt)
            << "\n";



//        output_l2_norm << vcl.getProcessingUnits() << "," << "rk4" << "," << std::accumulate(norm_l2_results.begin(), norm_l2_results.end(), 0.0) / static_cast<double>(norm_l2_results.size()) << "\n";


    }


    openfpm_finalize();

}
