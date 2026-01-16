

//! \cond [inclusion] \endcond
#include "Operators/Vector/vector_dist_operators.hpp"
#include "DCPSE/DCPSE_op/DCPSE_op.hpp"
#include <string>
//! \cond [inclusion] \endcond

constexpr int x = 0;
constexpr int y = 1;

typedef texp_v<double> state_type;

int main(int argc, char* argv[])
{



    // initialize the library
    openfpm_init(&argc,&argv);
    auto &v_cl=create_vcluster();



    //std::cout << "RK4"<< std::endl;

    std::ofstream output_file;
    output_file.open("rk4_error_results.csv",std::ios_base::app);

    size_t gdsz = std::stof(argv[1]);
    const size_t sz[3] = {gdsz,gdsz,gdsz};
    Box<3,double> box({0.0,0.0,0.0},{2.5,2.5,2.5});
    // Define periodicity of the grid
    size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};
    double spacing[3];
    spacing[0] = 2.5 / (sz[0]);
    spacing[1] = 2.5 / (sz[1]);
    spacing[2] = 2.5 / (sz[2]);
    double rCut = 2.9 * spacing[0];
    Ghost<3, double> ghost(rCut);

    // integration time
    //        double t=0.0,tf=0.5;
    double t=0,tf=std::stof(argv[2]);
    double init_dt = 1;

    typedef aggregate<VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>> Property_type;
    typedef vector_dist<3, double, Property_type> dist_vector_type;

    dist_vector_type Particles(0, box, bc, ghost);
    Particles.setPropNames({"Concentration","Initial","k1","k2","k3","k4"});
    // analytical solution
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


    auto verletList = Particles.template getVerlet<VL_NON_SYMMETRIC|VL_SKIP_REF_PART>(rCut);
    Laplacian<decltype(verletList)> Lap(Particles, verletList, 2, rCut, support_options::RADIUS);

    auto C = getV<0>(Particles);
    auto Init = getV<1>(Particles);
    auto k1 = getV<2>(Particles);
    auto k2 = getV<3>(Particles);
    auto k3 = getV<4>(Particles);
    auto k4 = getV<5>(Particles);
    // integrate for different numbers of time steps
    double dt = init_dt;

    //std::cout << std::endl;
    //std::cout << "--- RK4 Solving for dt = " << dt << " ..." << std::endl;
    size_t steps = 0;
    double K = 0.053;
    double F = 0.014;

    double d1 = 2*1e-4;
    double d2 = 1*1e-4;

    timer timer_integrate;
    timer_integrate.start();    
    while (t < tf) {
        C=Init;
        Particles.ghost_get<0>(SKIP_LABELLING);
        k1[x]=d1*Lap(C[x]) - C[x] * C[y] * C[y] + F - F * C[x];
        k1[y]=d2*Lap(C[y]) + C[x] * C[y] * C[y] - (F+K) * C[y];

        C=Init+0.5*dt*k1;
        Particles.ghost_get<0>(SKIP_LABELLING);
        k2[x]=d1*Lap(C[x]) - C[x] * C[y] * C[y] + F - F * C[x];
        k2[y]=d2*Lap(C[y]) + C[x] * C[y] * C[y] - (F+K) * C[y];
        
        C=Init+0.5*dt*k2;
        Particles.ghost_get<0>(SKIP_LABELLING);
        k3[x]=d1*Lap(C[x]) - C[x] * C[y] * C[y] + F - F * C[x];
        k3[y]=d2*Lap(C[y]) + C[x] * C[y] * C[y] - (F+K) * C[y];

        C=Init+dt*k3;
        Particles.ghost_get<0>(SKIP_LABELLING);
        k4[x]=d1*Lap(C[x]) - C[x] * C[y] * C[y] + F - F * C[x];
        k4[y]=d2*Lap(C[y]) + C[x] * C[y] * C[y] - (F+K) * C[y];
        Init = Init +  dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        //if(steps%100==0){
        //Particles.deleteGhost();
        //Particles.write_frame("Test",steps,t);
        //Particles.ghost_get<0>();
        //}
        t += dt;
        steps++;
    }
    timer_integrate.stop();
    Particles.ghost_get<0,1>();
    Particles.deleteGhost();
    Particles.write("Final");

    double tt=timer_integrate.getwct();
    v_cl.sum(tt);
    v_cl.execute();

    if(v_cl.rank()==0)
    {output_file << steps <<","<<v_cl.size()<<"," << tt/v_cl.size() << "\n";}
    // Particles.write("particles_dt=" + boost::lexical_cast<std::string>(dt));
    output_file.close();


    openfpm_finalize();

}
