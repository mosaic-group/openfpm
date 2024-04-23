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
#include "config.h"
#include "Vector/vector_dist_subset.hpp"
#include "DCPSE/DCPSE_op/DCPSE_surface_op.hpp"
#include "DCPSE/DCPSE_op/DCPSE_Solver.hpp"
#include "util/SphericalHarmonics.hpp"
#include <util/PathsAndFiles.hpp>

std::string p_output;

constexpr int DF=0;
constexpr int F=1;
constexpr int ANADF=2;
constexpr int ERR=3;
constexpr int NORMAL=4;
constexpr int PCOORD=5;
constexpr int ANAG=6;
constexpr int ERRG=7;
constexpr int POS=8;
constexpr int CURVT=9;
constexpr int PCURV=10;


typedef aggregate<double,double,double,double,Point<3,double>,Point<3,double>,double,double,Point<3,double>,double[3][3],double[3]> prop;
openfpm::vector<std::string> propNames{{"Gauss","Mean","AnaMean","ErrorMean","Normal","PCOORD","AnaGauss","ErrorGauss","Pos","CurvatureTensor","PrincipalCurvature"}};
typedef vector_dist_ws<3,double,prop> vector_type;

int main(int argc, char * argv[]) {
  openfpm_init(&argc,&argv);
  auto & v_cl = create_vcluster();
  timer tt;
  tt.start();
  size_t n;
  unsigned int ord;
  double spL;
  double grid_spacing_surf,sc;
  double rCut;

  // Get command line arguments
  if (argc != 4) {
    std::cout << "Error: The executable requires the following arguments: number_grid_points" << std::endl;
    return 0;
  }
  else {
    n = std::stoi(argv[1]);
    ord = std::stoi(argv[2]);
    sc = std::stoi(argv[3]);
    //spL=std::stof(argv[2]);
    //grid_spacing_surf=std::stof(argv[3]);
    //rCut=std::stof(argv[4]);
  }

  // Files and directories
  std::string cwd{boost::filesystem::current_path().string()};
  std::string output_dir{cwd + "/output_Ellipsoid"};
  create_directory_if_not_exist(output_dir);
  p_output = output_dir + "/particles";
  size_t n_sp=n;
  // Domain
  double boxP1{-1.5}, boxP2{1.5};
  double boxSize{boxP2 - boxP1};
  size_t sz[3] = {n,n,n};
  double grid_spacing{boxSize/(sz[0]-1)};
  grid_spacing_surf=grid_spacing*sc;
  rCut=2.9 * grid_spacing_surf;
  double a=1.0,b=0.8,c=0.75;

  Box<3,double> domain{{boxP1,boxP1,boxP1},{boxP2,boxP2,boxP2}};
  size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
  Ghost<3,double> ghost{rCut + grid_spacing/8.0};
  
  // particles
  vector_type particles{0,domain,bc,ghost};
  particles.setPropNames(propNames);

  // 1. particles on the Ellipsoidal surface //phi=u theta=v
  double theta{0.0},phi{0.0};
  double dtheta{M_PI/double(n-1)};
  double dphi{2*M_PI/double(2.0*n-1)};
  double Golden_angle=M_PI * (3.0 - sqrt(5.0));
  if (v_cl.rank() == 0) {
    for(int i=0;i<n;i++)
      {for(int j=0;j<2*n-1;j++)
        {
        particles.add();
        double xp=a*std::sin(theta)*std::cos(phi);
        double yp=b*std::sin(theta)*std::sin(phi);
        double zp=c*std::cos(theta);
        particles.getLastPos()[0] = xp;
        particles.getLastPos()[1] = yp;
        particles.getLastPos()[2] = zp;
        particles.getLastProp<POS>()[0] = xp;
        particles.getLastProp<POS>()[1] = yp;
        particles.getLastProp<POS>()[2] = zp;
        double rm=sqrt(4*(xp*xp/(a*a*a*a)+yp*yp/(b*b*b*b)+zp*zp/(c*c*c*c)));
        particles.getLastProp<F>() = 0;
        particles.getLastProp<NORMAL>()[0] = 2*xp/(a*a*rm);
        particles.getLastProp<NORMAL>()[1] = 2*yp/(b*b*rm);
        particles.getLastProp<NORMAL>()[2] = 2*zp/(c*c*rm);
        particles.getLastProp<PCOORD>()[0] = 0;
        particles.getLastProp<PCOORD>()[1] = theta;
        particles.getLastProp<PCOORD>()[2] = phi;
        particles.getLastProp<ANADF>() = (a*b*c)*(3*(a*a+b*b)+2*c*c+(a*a+b*b-2*c*c)*std::cos(2*theta)-2*(a*a-b*b)*std::cos(2*phi)*std::sin(theta)*std::sin(theta))/(8*openfpm::math::intpowlog(sqrt((a*a*b*b*std::cos(theta)*std::cos(theta)+c*c*(b*b*std::cos(phi)*std::cos(phi)+a*a*std::sin(phi)*std::sin(phi))*std::sin(theta)*std::sin(theta))),3));
        particles.getLastProp<ANAG>() = (a*a*b*b*c*c)/(openfpm::math::intpowlog((a*a*b*b*std::cos(theta)*std::cos(theta)+c*c*(b*b*std::cos(phi)*std::cos(phi)+a*a*std::sin(phi)*std::sin(phi))*std::sin(theta)*std::sin(theta)),2));
        particles.getLastSubset(0);
        phi += dphi;
        if(i==0 || i==n-1)
        { break;}
        }
        theta+=dtheta;
        }
    std::cout << "n: " << n*2*n << " - grid spacing: " << grid_spacing << " - rCut: " << rCut << "Nspacing" << grid_spacing_surf<< " - dTheta: "<<dtheta<<" - dPhi: "<<dphi <<std::endl;
  }

  particles.map();
  particles.ghost_get<F>();
  particles.deleteGhost();
  particles.write_frame("p_output",BINARY);
  vector_dist_subset<3,double,prop> particles_bulk(particles,0);
  vector_dist_subset<3,double,prop> particles_boundary(particles,1);
  auto &bulkIds=particles_bulk.getIds();
  auto &bdrIds=particles_boundary.getIds();
  
  auto DErr=getV<ERR>(particles);
  auto Af=getV<ANADF>(particles);
  auto Df=getV<DF>(particles);
  auto f=getV<F>(particles);
  auto N=getV<NORMAL>(particles);
  auto Pos=getV<POS>(particles);
  auto CurvTensor=getV<CURVT>(particles);
  particles.ghost_get<F,NORMAL,POS>();
  auto verletList = particles.getVerletWithoutRefP(rCut);
  unsigned int Ncount= (rCut/grid_spacing_surf);
  SurfaceDerivative_x<NORMAL> Sdx{particles,ord,verletList,grid_spacing_surf,Ncount};
  SurfaceDerivative_y<NORMAL> Sdy{particles,ord,verletList,grid_spacing_surf,Ncount};
  SurfaceDerivative_z<NORMAL> Sdz{particles,ord,verletList,grid_spacing_surf,Ncount};

  auto Na=Sdx(N[0]),Nb=Sdy(N[0]),Nc=Sdz(N[0]),
       Nd=Sdx(N[1]),Ne=Sdy(N[1]),Nf=Sdz(N[1]),
       Ng=Sdx(N[2]),Nh=Sdy(N[2]),Ni=Sdz(N[2]);
  //SurfaceDerivative_xx<NORMAL> Sdxx{particles,2,rCut,grid_spacing_surf};
  //SurfaceDerivative_yy<NORMAL> Sdyy{particles,2,rCut,grid_spacing_surf};
  //SurfaceDerivative_zz<NORMAL> Sdzz{particles,2,rCut,grid_spacing_surf};
  f=(Na+Ne+Ni)/2.0;  //2H=-Div(n_hat) //Note that sign of the normal matters. in our case it is negative.
  //f=-0.5*((Sdxx(Pos[0])+Sdyy(Pos[0])+Sdzz(Pos[0]))*N[0]+(Sdxx(Pos[1])+Sdyy(Pos[1])+Sdzz(Pos[1]))*N[1]+(Sdxx(Pos[2])+Sdyy(Pos[2])+Sdzz(Pos[2]))*N[2]);
  //Df=(Na*Ne*Ni - Na*Nf*Nh - Nb*Nd*Ni + Nb*Nf*Ng + Nc*Nd*Nh - Nc*Ne*Ng); This doesnt work because there is a 0 eigen value

  CurvTensor[0][0]=Na;
  CurvTensor[0][1]=Nb;
  CurvTensor[0][2]=Nc;
  CurvTensor[1][0]=Nd;
  CurvTensor[1][1]=Ne;
  CurvTensor[1][2]=Nf;
  CurvTensor[2][0]=Ng;
  CurvTensor[2][1]=Nh;
  CurvTensor[2][2]=Ni;

  auto it=particles.getDomainIterator();
  while(it.isNext())
  {
    auto p=it.get();
    typedef EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    MatrixType CT(3, 3);
    for(int i=0;i<3;i++)
    {
      for(int j=0;j<3;j++)
      {CT(i, j)= particles.getProp<CURVT>(p)[i][j];}
    }
    //Eigen::EigenSolver<MatrixType> es(CT);
    double a=CT.eigenvalues()[0].real(),b=CT.eigenvalues()[1].real(),c=CT.eigenvalues()[2].real();
    if (a > c) std::swap(a, c);
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    //std::cout<<"\nEigs: "<<CT.eigenvalues();
    particles.getProp<PCURV>(p)[0]=c;
    particles.getProp<PCURV>(p)[1]=b;
    particles.getProp<PCURV>(p)[2]=a;
    particles.getProp<DF>(p)=b*c;
    ++it;
  }

  particles.deleteGhost();
  particles.write_frame(p_output,0,BINARY);
  /*DCPSE_scheme<equations3d1, decltype(particles)> Solver(particles);
  auto Eig = Sdxx(f)+Sdyy(f)+Sdzz(f)+K*(K+1)*f;
  Solver.impose(Eig, bulkIds, 0);
  Solver.impose(f, bdrIds, Af);
  Solver.solve(f);*/

  // Error
  double MaxError=0,L2=0,MaxErrorG=0,L2G=0;
  for (int j = 0; j < bulkIds.size(); j++) {
      auto p = bulkIds.get<0>(j);
      particles.getProp<ERR>(p) = fabs(particles.getProp<F>(p)-particles.getProp<ANADF>(p));
      particles.getProp<ERRG>(p) = fabs(particles.getProp<DF>(p)-particles.getProp<ANAG>(p));
      if (particles.getProp<ERR>(p)>MaxError){
              MaxError=particles.getProp<ERR>(p);
          }
      if (particles.getProp<ERRG>(p)>MaxErrorG){
              MaxErrorG=particles.getProp<ERRG>(p);
          }
          L2+=particles.getProp<ERR>(p)*particles.getProp<ERR>(p);
          L2G+=particles.getProp<ERRG>(p)*particles.getProp<ERRG>(p);
  }
  v_cl.sum(L2);
  v_cl.sum(L2G);
  v_cl.max(MaxError);
  v_cl.max(MaxErrorG);
  v_cl.execute();
  L2=sqrt(L2);
  L2G=sqrt(L2G);
  std::cout.precision(16);
  if (v_cl.rank()==0){
      std::cout<<"Mean Curvature L2:"<<L2<<std::endl;
      std::cout<<"Mean Curvature L_inf:"<<MaxError<<std::endl;
      std::cout<<"Gauss Curvature L2:"<<L2G<<std::endl;
      std::cout<<"Gauss Curvature L_inf:"<<MaxErrorG<<std::endl;
  }
  particles.deleteGhost();
  particles.write(p_output,BINARY);
  //Sdxx.deallocate(particles);
  //Sdyy.deallocate(particles);
  //Sdzz.deallocate(particles);

  tt.stop();
  if (v_cl.rank() == 0)
    std::cout << "Simulation took: " << tt.getcputime() << " s (CPU) - " << tt.getwct() << " s (wall)\n";
  
  openfpm_finalize();
  return 0;

}
