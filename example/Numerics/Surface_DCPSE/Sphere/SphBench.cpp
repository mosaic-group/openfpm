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

//#define SE_CLASS1

#include "Vector/vector_dist_subset.hpp"
#include "DCPSE/DCPSE_op/DCPSE_surface_op.hpp"
//#include "DCPSE/DCPSE_op/DCPSE_Solver.hpp"
#include "util/SphericalHarmonics.hpp"
#include <util/PathsAndFiles.hpp>

std::string p_output;

constexpr int DF=0;
constexpr int F=1;
constexpr int ANADF=2;
constexpr int ERR=3;
constexpr int NORMAL=4;
constexpr int PCOORD=5;

typedef aggregate<double,double,double,double,Point<3,double>,Point<3,double>> prop;
openfpm::vector<std::string> propNames{{"Df","f","AnaDf","error","Normal","PCOORD"}};
typedef vector_dist_ws<3,double,prop> vector_type;

int main(int argc, char * argv[]) {
  openfpm_init(&argc,&argv);
  auto & v_cl = create_vcluster();
  timer tt,tt2,tt3,tt4;
  tt4.start();
  size_t n;
  double spL;
  double grid_spacing_surf;
  double rCut;
  constexpr int K = 4;

  // Get command line arguments
  if (argc > 3) {
    if(v_cl.rank()==0)
    std::cout << "Warning: The executable requires the following arguments: number_grid_points orderOfConv" << std::endl;
  }
  
  n = std::stoi(argv[1]);
  unsigned int ord=std::stoi(argv[2]);
    //spL=std::stof(argv[2]);
    //grid_spacing_surf=std::stof(argv[3]);
    //rCut=std::stof(argv[4]);
   
  p_output = "particles";
  size_t n_sp=n;
  // Domain
  double boxP1{-1.5}, boxP2{1.5};
  double boxSize{boxP2 - boxP1};
  size_t sz[3] = {n,n,n};
  double grid_spacing{0.8/((std::pow(sz[0],1/3.0)-1))};
  double Golden_angle=M_PI * (3.0 - sqrt(5.0));
  grid_spacing_surf=grid_spacing;
  rCut=2.9 * grid_spacing_surf;

  Box<3,double> domain{{boxP1,boxP1,boxP1},{boxP2,boxP2,boxP2}};
  size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
  Ghost<3,double> ghost{rCut + grid_spacing/8.0};
  
  // particles
  vector_type particles{0,domain,bc,ghost};
  particles.setPropNames(propNames);

  Point<3,double> coord;
  // 1. particles on the Circular surface
  double theta{0.0};
  if (v_cl.rank() == 0) {
    for(int i=1;i<n_sp;i++)
        {
            double y = 1.0 - (i /double(n_sp - 1.0)) * 2.0;
            double radius = sqrt(1 - y * y);
            double Golden_theta = Golden_angle * i;
            double x = cos(Golden_theta) * radius;
            double z = sin(Golden_theta) * radius;
            particles.add();
            particles.getLastPos()[0] = x;
            particles.getLastPos()[1] = y;
            particles.getLastPos()[2] = z;
            double rm=sqrt(x*x+y*y+z*z);
            particles.getLastProp<NORMAL>()[0] = x/rm;
            particles.getLastProp<NORMAL>()[1] = y/rm;
            particles.getLastProp<NORMAL>()[2] = z/rm;
            particles.getLastProp<PCOORD>()[0] = 1.0 ;
            particles.getLastProp<PCOORD>()[1] = std::atan2(sqrt(x*x+y*y),z);
            particles.getLastProp<PCOORD>()[2] = std::atan2(y,x);
        }
        particles.getLastSubset(1);
    std::cout << "n: " << n <<" - rCut: " << rCut << " nSpacing" << grid_spacing_surf<<std::endl;
  }

  particles.map();
  particles.ghost_get<F>();

  vector_dist_subset<3,double,prop> particles_bulk(particles,0);
  vector_dist_subset<3,double,prop> particles_boundary(particles,1);
  auto &bulkIds=particles_bulk.getIds();
  auto &bdrIds=particles_boundary.getIds();
  std::unordered_map<const lm,double,key_hash,key_equal> Alm;
  //Setting max mode l_max
  //Setting amplitudes to 1
  for(int l=0;l<=K;l++){
      for(int m=-l;m<=l;m++){
          Alm[std::make_tuple(l,m)]=0;
      }
  }
  //2Alm[std::make_tuple(1,0)]=1;
  Alm[std::make_tuple(K,0)]=1;
  //Alm[std::make_tuple(3,2)]=1;

  spL=K;
  auto it2 = particles.getDomainIterator();
  while (it2.isNext()) {
      auto p = it2.get();
      Point<3, double> xP = particles.getProp<PCOORD>(p);
      particles.getProp<ANADF>(p)=openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      particles.getProp<DF>(p)=-(K)*(K+1)*openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      ++it2;
  }
  particles.ghost_get<ANADF>();
  auto DErr=getV<ERR>(particles);
  auto Af=getV<ANADF>(particles);
  auto Df=getV<DF>(particles);
  auto f=getV<F>(particles);
  tt.start();
  tt3.start();
   auto verletList = particles.getVerletWithoutRefP(rCut);
  unsigned int nCount=(rCut/grid_spacing_surf);
  SurfaceDerivative_xx<NORMAL> Sdxx{particles,ord,verletList,grid_spacing_surf,nCount};
  SurfaceDerivative_yy<NORMAL> Sdyy{particles,ord,verletList,grid_spacing_surf,nCount};
  SurfaceDerivative_zz<NORMAL> Sdzz{particles,ord,verletList,grid_spacing_surf,nCount};
  tt3.stop();
  tt2.start();
  particles.ghost_get<F>();
  f=(Sdxx(Af)+Sdyy(Af)+Sdzz(Af));
  tt2.stop();
  tt.stop();
  particles.deleteGhost();
  particles.write_frame(p_output,0,BINARY);
  // Error
  double MaxError=0,L2=0;
  for (int j = 0; j < bulkIds.size(); j++) {
      auto p = bulkIds.get<0>(j);
      particles.getProp<ERR>(p) = fabs(particles.getProp<F>(p)-particles.getProp<DF>(p));
      if (particles.getProp<ERR>(p)>MaxError){
              MaxError=particles.getProp<ERR>(p);
          }
          L2+=particles.getProp<ERR>(p)*particles.getProp<ERR>(p);
  }
  double tcpu=tt.getcputime(),twall=tt.getwct(),
         Dtcpu=tt2.getcputime(),Dtwall=tt2.getwct(),
         DCtcpu=tt3.getcputime(),DCtwall=tt3.getwct();
  v_cl.sum(DCtcpu);
  v_cl.sum(DCtwall);
  v_cl.sum(Dtcpu);
  v_cl.sum(Dtwall);
  v_cl.sum(tcpu);
  v_cl.sum(twall);
  v_cl.sum(L2);
  v_cl.max(MaxError);
  v_cl.execute();
  L2=sqrt(L2);
  std::cout.precision(16);
  if (v_cl.rank()==0){
      std::cout<<"L2:"<<L2/double(sqrt(n_sp))<<std::endl;
      std::cout<<"L_inf:"<<MaxError<<std::endl;
  }

  particles.deleteGhost();
  particles.write(p_output,BINARY);
  Sdxx.deallocate(particles);
  Sdyy.deallocate(particles);
  tt4.stop();
  double Ftcpu=tt4.getcputime(),Ftwall=tt.getwct();
  v_cl.sum(Ftcpu);
  v_cl.sum(Ftwall);
  v_cl.execute();
  if (v_cl.rank() == 0){
    std::cout << "DCPSE Construction took: " << DCtcpu/v_cl.size() << " s (CPU) - " << DCtwall/v_cl.size() << " s (wall)\n";
    std::cout << "DCPSE Evaluation took: " << Dtcpu/v_cl.size() << " s (CPU) - " << Dtwall/v_cl.size() << " s (wall)\n";
    std::cout << "Total DCPSE time: " << tcpu/v_cl.size() << " s (CPU) - " << twall/v_cl.size() << " s (wall)\n";
    std::cout << "Entire Simulation took: " << Ftcpu/v_cl.size() << " s (CPU) - " << Ftwall/v_cl.size() << " s (wall)\n";
  }
  openfpm_finalize();
  return 0;

}
