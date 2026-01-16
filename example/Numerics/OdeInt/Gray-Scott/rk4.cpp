

//! \cond [inclusion] \endcond
#include <stddef.h>
#include "Vector/vector_dist.hpp"
#include <iostream>
#include <fstream>
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "Decomposition/Distribution/SpaceDistribution.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"
//! \cond [inclusion] \endcond



typedef texp_v<double> state_type;

void Exponential2( const state_type &x , state_type &dxdt , const double t )
{
    dxdt = x;
}

void sigmoid2( const state_type &x , state_type &dxdt , const double t )
{
    // g'(t) = g(t) * (1-sigma)
    dxdt = x * (1- (1/(1+exp(-t))) );
}

int main(int argc, char* argv[])
{



    // initialize the library
    openfpm_init(&argc,&argv);
    auto &v_cl=create_vcluster();



    //std::cout << "RK4"<< std::endl;

    std::ofstream output_file;
    output_file.open("rk4_error_results.csv",std::ios_base::app);

    size_t edgeSemiSize = std::stof(argv[1]);
    const size_t sz[2] = {edgeSemiSize,edgeSemiSize };
    Box<2, double> box({ 0, 0 }, { 1.0, 1.0 });
    size_t bc[2] = { NON_PERIODIC, NON_PERIODIC };
    double spacing[2];
    spacing[0] = 1.0 / (sz[0] - 1);
    spacing[1] = 1.0 / (sz[1] - 1);
    double rCut = 3.9 * spacing[0];
    Ghost<2, double> ghost(rCut);

    // integration time
    //        double t=0.0,tf=0.5;
    double t=-5.0,tf=5.0;
    double init_dt = 0.01;

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

    state_type x0;

    // integrate for different numbers of time steps

    double dt = init_dt;

    //std::cout << std::endl;
    //std::cout << "--- RK4 Solving for dt = " << dt << " ..." << std::endl;

    x0 = Init;

    size_t steps = 0;


    timer timer_integrate;
    timer_integrate.start();

    //            steps = integrate_adaptive( dopri5_controlled_stepper_type() , Exponential2 , x0 , t , tf , init_dt);
    // steps = boost::numeric::odeint::integrate_const(rk4_stepper_type(), sigmoid2, x0, t, tf, dt);


    void (*derivative)(const state_type& , state_type& , const double);
    derivative = &Exponential2;

    while (t < tf - dt) {
        state_type k1;
        state_type k2;
        state_type k3;
        state_type k4;

        derivative(x0, k1, t);

        derivative(x0 + 0.5 * dt * k1, k2, t + 0.5 * dt);

        derivative(x0 + 0.5 * dt * k2, k3, t + 0.5 * dt);

        derivative(x0 + dt * k3, k4, t + dt);

        x0 = x0 +  dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);

        t += dt;
        steps++;
    }

    timer_integrate.stop();

    OdeSol = x0;

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

    auto & vcl = create_vcluster();
    if (vcl.getProcessUnitID() == 0)
    {
    std::cout << "Steps: " << steps << std::endl;
    std::cout << "Time: " << timer_integrate.getwct() << std::endl;
    std::cout << "Norm inf: " << norm_inf << std::endl;
    std::cout << "Norm L2: " << norm_l2 << std::endl;
    }
    double tt=timer_integrate.getwct();
    v_cl.sum(tt);
    v_cl.execute();

    if(v_cl.rank()==0)
    {output_file << steps <<","<<v_cl.size()<<"," << tt/v_cl.size() << "," << norm_inf << "," << norm_l2 << "\n";}
    // Particles.write("particles_dt=" + boost::lexical_cast<std::string>(dt));
    output_file.close();


    openfpm_finalize();

}