// ------------------------------------------------------------------------------------------
//    Copyright (C) 2021 ---- absingh
//
//    This file is part of the Surface DCPSE Paper.
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ------------------------------------------------------------------------------------------

#include "Vector/vector_dist_subset.hpp"
#include "DCPSE/DCPSE_op/DCPSE_surface_op.hpp"
#include <util/PathsAndFiles.hpp>
#include "OdeIntegrators/OdeIntegrators.hpp"

void *PointerGlobal, *Pointer2Bulk,*Pointer2Boundary;

typedef aggregate<double,double,Point<3,double>,double,double>  prop;
typedef vector_dist_ws<3,double,prop> vector_type;
typedef vector_dist_subset<3,double,prop> subset_type;

constexpr int CONC=0;
constexpr int DCONC=1;
constexpr int NORMAL=2;
constexpr int COLD=3;
openfpm::vector<std::string> propNames{{"Concentration","DConc","Normal","Conc_old","Error"}};

double dt,tf,max_steady_tol;
int wr_at;
constexpr int dim=3;

template<int dim>
struct bump {
  const Point<2,double>  & center;
  const double alpha;
  const double radius;
  const double threshold;
};

template <int dim>
double f( Point<3,double>  coord,
bump<dim> & surf) {

double arg, arg2;

arg = 0.0;
for (int d = 0; d < dim-1; ++d)
arg += (coord[d]-surf.center[d])*(coord[d]-surf.center[d]);
arg2 = arg/(surf.radius*surf.radius);

if (arg2 < (1 - surf.threshold)*(1 - surf.threshold))
return surf.alpha * std::exp(-1.0/(1.0-arg2));
else
return 0.0;
}

template<int dim>
Point<3,double> init_normal(Point<3,double> coord,
                   bump<dim> & surf) {

  Point<3,double> normal;
  double x, y, dist, arg, prefactor, norm_grad;
  const double r2{surf.radius*surf.radius};
  
  dist = 0.0;
  for (int d = 0; d < dim-1; ++d) 
    dist += (coord[d]-surf.center[d])*(coord[d]-surf.center[d]);

  if (std::fabs(dist - r2) < 1e-10)
    arg = r2 / 1e-10;
  else
    arg = r2 / (dist - r2);
  
  x = coord[0]-surf.center[0];
  y = coord[1]-surf.center[1];
  prefactor = - 2.0 * f<dim>(coord,surf) * arg*arg / r2;
  
  norm_grad = std::sqrt(4.0 * f<dim>(coord,surf)*f<dim>(coord,surf) * arg*arg*arg*arg * dist / (r2*r2) + 1.0);

  normal[0] = - prefactor * x / norm_grad;
  normal[1] = - prefactor * y / norm_grad;
  normal[2] = 1.0 / norm_grad;

  double mag=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  normal=normal/mag;
  
  return normal;
}

Point<3,double> SurfaceNormal(double &x,double &y,double &z,double alpha){
    Point<3,double> Normal={0,0,1};
    double dist=sqrt((x+0.5)*(x+0.5)+y*y);
    if(dist>(1.0-0.025))
    {   
        return Normal;
    }
    else{
        dist=dist/0.25;
        double norm_grad = -alpha*2*dist*exp(-1./(1.-dist*dist)*1./((1.-dist*dist)*(1.-dist*dist)));
        Normal[0]=x;
        Normal[1]=y;
        Normal[2]=4.0*z*norm_grad*(0.5/dist)*(2*(x+0.5)+2*y);
        return Normal/Normal.distance(Normal);
    }
    }

template<typename DXX,typename DYY,typename DZZ>
struct RHSFunctor
{
    //Intializing the operators
    DXX &SDxx;
    DYY &SDyy;
    DZZ &SDzz;
    //Constructor
    RHSFunctor(DXX &SDxx,DYY &SDyy,DZZ &SDzz):SDxx(SDxx),SDyy(SDyy),SDzz(SDzz)
    {}

    void operator()( const state_type_1d_ofp &X , state_type_1d_ofp &dxdt , const double t ) const
    {
        //Casting the pointers to OpenFPM vector distributions
        vector_type &Particles= *(vector_type *) PointerGlobal;
        subset_type &Particles_bulk= *(subset_type *) Pointer2Bulk;

        //Aliasing the properties.
        auto C = getV<CONC>(Particles);
        auto C_bulk = getV<CONC>(Particles_bulk);
        auto dC = getV<DCONC>(Particles);
        auto dC_bulk = getV<DCONC>(Particles_bulk);
        C_bulk=X.data.get<0>();
        Particles.ghost_get<CONC>();

        // We do the RHS computations for the Laplacian (Updating bulk only).
        dC_bulk = (SDxx(C)+SDyy(C)+SDzz(C));

        //We copy back to the dxdt state_type for Odeint
        dxdt.data.get<0>()=dC;
    }
};

struct ObserverFunctor {
    int ctr;
    int wctr; 
    double t_old;
    //Constructor
    ObserverFunctor(){
        //a counter for counting the np. of steps
        ctr = 0;
        wctr= 0;
        //Starting with t=0, we compute the step size take by t-t_old. So for the first observed step size is what we provide. Which will be 0-(-dt)=dt.
        t_old = -dt;
    }
    void operator()(state_type_1d_ofp &X, double t) {
        //Casting the pointers to OpenFPM vector distributions
        auto & v_cl = create_vcluster();
        vector_type &Particles = *(vector_type *) PointerGlobal;
        subset_type &Particles_bulk = *(subset_type *) Pointer2Bulk;
        auto &bulk=Particles_bulk.getIds();
        auto Concentration = getV<CONC>(Particles);
        auto Concentration_old = getV<COLD>(Particles);
        auto Concentration_bulk = getV<CONC>(Particles_bulk);
        Concentration_bulk = X.data.get<0>();
        if (v_cl.rank() == 0) {
                        std::cout << "Time step " << ctr << " : " << t << " over." <<"dt is set to: "<<(t-t_old)<< std::endl;
                        std::cout << "----------------------------------------------------------" << std::endl;
                    }
        if(ctr%wr_at==0){
          Particles.deleteGhost();
          Particles.write_frame("ERDiff", wctr,t,BINARY);
          Particles.ghost_get<0>();
          wctr++;
        }
        double MaxRateOfChange=0;
        for (int j = 0; j < bulk.size(); j++) {
            auto p = bulk.get<0>(j);
            if(Particles.getProp<CONC>(p)<0){
                Particles.getProp<CONC>(p)=0;
            }
            if(fabs(Particles.getProp<CONC>(p)-Particles.getProp<COLD>(p))>MaxRateOfChange)
            {
                MaxRateOfChange=fabs(Particles.getProp<CONC>(p)-Particles.getProp<COLD>(p));
            }
        }
        v_cl.max(MaxRateOfChange);
        v_cl.execute();
        if(v_cl.rank()==0)
        {std::cout<<"MaxRateOfChange: "<<MaxRateOfChange<<std::endl;
        }
        if(MaxRateOfChange<max_steady_tol && ctr>5)
        {
            std::cout<<"Steady State Reached."<<std::endl;
            openfpm_finalize();
            exit(0);
        }
        Concentration_old=Concentration;
        Particles.ghost_get<0>();
        ctr++;
        t_old=t;
        X.data.get<0>()=Concentration;
    }
};

int main(int argc, char * argv[]) {
  openfpm_init(&argc,&argv);
  auto & v_cl = create_vcluster();
  timer tt;
  tt.start();
  double grid_spacing_surf;
  double rCut,SCF,alph;
  bool DCPSE_LOAD;

  // Get command line arguments
  std::ifstream PCfile;
  if (argc < 9) {
    std::cout << "Error: Not exact no. of args." << std::endl;
    return 0;
  }
  else {
    grid_spacing_surf=0.03125;
    PCfile.open(argv[1]);
    tf=std::stof(argv[2]);
    dt=std::stof(argv[3]);
    wr_at=std::stoi(argv[4]);
    max_steady_tol=std::stof(argv[5]);
    SCF=std::stof(argv[6]);
    alph=std::stof(argv[7]);
    DCPSE_LOAD=std::stof(argv[8]);
  }

  rCut=3.1*grid_spacing_surf;

  Box<3,double> domain{{-2.5,-2.5,-2.5},{2.5,2.5,2.5}};
  size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
  Ghost<3,double> ghost{rCut + grid_spacing_surf/8.0};
  struct bump<dim> surf = {.center= Point<dim-1,double>{-0.5,0.0}, .alpha=alph, .radius=0.25, .threshold=0.025};

    // particles
  vector_type particles{0,domain,bc,ghost};
  particles.setPropNames(propNames);
  double uconc,Nx,Ny,Nz,Px,Py,Pz;
  if(v_cl.rank()==0){
    int pctr=0;
    while ( PCfile >>uconc>> Px >> Py >> Pz)
      {
          particles.add();
          particles.getLastPos()[0]=Px;
          particles.getLastPos()[1]=Py;
          particles.getLastPos()[2]=Pz;
          particles.getLastProp<CONC>()=uconc;
          particles.getLastProp<COLD>()=sin(2*Px)*sin(2*Py);
          particles.getLastProp<DCONC>()=-8*sin(2*Px)*sin(2*Py);
          particles.getLastProp<NORMAL>()=init_normal<dim>(particles.getLastPos(),surf);

          particles.getLastSubset(0);
          if(Px<-2+1e-3){
            particles.getLastSubset(1);
          }
          else if(Px>2-1e-3){
            particles.getLastSubset(1);
          }
          else if(Py<-2+1e-3){
            particles.getLastSubset(1);
          }
          else if(Py>2-1e-3){
            particles.getLastSubset(1);
          }

          ++pctr;
      }
    std::cout << "n: " << pctr << " - rCut: " << rCut << " - nSpacing: " << grid_spacing_surf<<" - Final Time: "<<tf<<" - Initial DT: "<<dt<<std::endl;
  }

  particles.map();
  particles.ghost_get<CONC>();
  particles.ghost_get_subset();
  particles.write("Init",BINARY);
  vector_dist_subset<3,double,prop> particles_bulk(particles,0);
  vector_dist_subset<3,double,prop> particles_boundary(particles,1);
  auto &bulk=particles_boundary.getIds();
  auto C=getV<CONC>(particles);
  auto DC=getV<DCONC>(particles);
  auto C_old=getV<COLD>(particles);


  //DC=0;
  timer ttt;
        ttt.start();
        create_directory_if_not_exist("DCPSE");
        support_options opt;
        if(DCPSE_LOAD){
            opt=support_options::LOAD;            
        }
        else
        {
            opt=support_options::ADAPTIVE;
        }
        auto verletList = particles.template getVerlet<VL_NON_SYMMETRIC|VL_SKIP_REF_PART>(rCut);
        SurfaceDerivative_xx<NORMAL,decltype(verletList)> Sdxx{particles,verletList,2,rCut,SCF,rCut/SCF,opt};
        SurfaceDerivative_yy<NORMAL,decltype(verletList)> Sdyy{particles,verletList,2,rCut,SCF,rCut/SCF,opt};
        SurfaceDerivative_zz<NORMAL,decltype(verletList)> Sdzz{particles,verletList,2,rCut,SCF,rCut/SCF,opt};
        if(DCPSE_LOAD){
            Sdxx.load(particles,"DCPSE/Dxx");
            Sdyy.load(particles,"DCPSE/Dyy");
            Sdzz.load(particles,"DCPSE/Dzz");
        }
        else{
            Sdxx.save(particles,"DCPSE/Dxx");
            Sdyy.save(particles,"DCPSE/Dyy");
            Sdzz.save(particles,"DCPSE/Dzz");
        }
        /*DC=0;
        C=0;
        C_old=0;
        auto itNNN=particles.getDomainIterator();
        while(itNNN.isNext()){
            auto p=itNNN.get().getKey();
            Point<3,double> xp=particles.getPos(p);
            if(xp.distance(0)<1e-2)
            {
            Sdxx.DrawKernel<CONC,vector_type>(p);
            Sdyy.DrawKernel<COLD,vector_type>(p);
            Sdzz.DrawKernel<DCONC,vector_type>(p);
            particles.write_frame("Kernel",p);
            DC=0;
            C=0;
            C_old=0;
            }
            ++itNNN;
        }*/
        ttt.stop();
        if(v_cl.rank()==0)
            std::cout<<"DCPSE Time: "<<tt.getwct()<<" seconds."<<std::endl;
  PointerGlobal = (void *) &particles;
  Pointer2Bulk = (void *) &particles_bulk;
  Pointer2Boundary = (void *) &particles_boundary;
  RHSFunctor<SurfaceDerivative_xx<NORMAL,decltype(verletList)>, SurfaceDerivative_yy<NORMAL,decltype(verletList)>,SurfaceDerivative_zz<NORMAL,decltype(verletList)>> System(Sdxx,Sdyy,Sdzz);
  ObserverFunctor Obs;
  state_type_1d_ofp X;
  X.data.get<0>() = C;

 particles.ghost_get<COLD>();
 //C=Sdxx(C_old)+Sdyy(C_old)+Sdzz(C_old);

  double worst1 = 0.0;

    for(int j=0;j<bulk.size();j++)
    {   auto p=bulk.get<0>(j);
        if (fabs(particles.getProp<DCONC>(p) - particles.getProp<CONC>(p)) >= worst1) {
            worst1 = fabs(particles.getProp<DCONC>(p) - particles.getProp<CONC>(p));
        }
        particles.getProp<4>(p) = fabs(particles.getProp<DCONC>(p) - particles.getProp<CONC>(p));

    }
    std::cout << "Maximum Analytic Error: " << worst1 << std::endl;
    DC=0;
  //boost::numeric::odeint::adaptive_adams_bashforth_moulton<2, state_type_1d_ofp,double,state_type_1d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp > abmA;
  //size_t steps=boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(max_steady_tol*0.1,max_steady_tol*0.1,abmA),System,X,0.0,tf,dt,Obs);
    boost::numeric::odeint::euler< state_type_1d_ofp,double,state_type_1d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp > Euler;
    size_t steps=boost::numeric::odeint::integrate_const(Euler,System,X,0.0,tf,dt,Obs);

  particles.deleteGhost();
  particles.write("Final",BINARY);
  Sdxx.deallocate(particles);
  Sdyy.deallocate(particles);
  Sdzz.deallocate(particles);

  tt.stop();
  if (v_cl.rank() == 0)
    std::cout << "Simulation took: " << tt.getcputime() << " s (CPU) - " << tt.getwct() << " s (wall)\n";
  
  openfpm_finalize();
  return 0;

}