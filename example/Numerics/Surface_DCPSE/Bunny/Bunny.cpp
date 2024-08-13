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

typedef aggregate<double,double,double,double,Point<3,double>,Point<3,double>> prop;
openfpm::vector<std::string> propNames{{"Df","f","AnaDf","error","Normal","PCOORD"}};
typedef vector_dist_ws<3,double,prop> vector_type;

int main(int argc, char * argv[]) {
  openfpm_init(&argc,&argv);
  auto & v_cl = create_vcluster();
  timer tt;
  tt.start();
  double grid_spacing_surf;
  double rCut;
  double EIGENVAL;

  // Get command line arguments
  if (argc < 3) {
    std::cout << "Error: The executable requires the following arguments: number_grid_points" << std::endl;
    return 0;
  }
  else {
    //n = std::stoi(argv[1]);
    //spL=std::stof(argv[2]);
    grid_spacing_surf=std::stof(argv[1]);
    EIGENVAL=std::stof(argv[2]);
    //rCut=std::stof(argv[4]);
  }

  // Files and directories
  std::string cwd{boost::filesystem::current_path().string()};
  std::string output_dir{cwd + "/output_Bunny"};
  create_directory_if_not_exist(output_dir);
  p_output = output_dir + "/particles";
  double boxP1{-1.5}, boxP2{1.5};
  double boxSize{boxP2 - boxP1};
  rCut=2.5 * grid_spacing_surf;

  Box<3,double> domain{{boxP1,boxP1,boxP1},{boxP2,boxP2,boxP2}};
  size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
  Ghost<3,double> ghost{rCut + grid_spacing_surf/8.0};
  
    // particles
  vector_type particles{0,domain,bc,ghost};
  particles.setPropNames(propNames);
  double Nx,Ny,Nz,Px,Py,Pz;
  std::ifstream PCfile;
  PCfile.open("data.csv");
  if(v_cl.rank()==0){
    int pctr=0;
    while ( PCfile >> Nx >> Ny >> Nz >> Px >> Py >> Pz )
      {
          particles.add();
          particles.getLastPos()[0]=Px;
          particles.getLastPos()[1]=Py;
          particles.getLastPos()[2]=Pz;
          particles.getLastProp<NORMAL>()[0]=Nx;
          particles.getLastProp<NORMAL>()[1]=Ny;
          particles.getLastProp<NORMAL>()[2]=Nz;
          ++pctr;
      }
    particles.getLastSubset(1);
    std::cout << "n: " << pctr << " - rCut: " << rCut << "nSpacing" << grid_spacing_surf<<std::endl;
  }

  particles.map();
  particles.ghost_get<F>();
  particles.ghost_get_subset();

  vector_dist_subset<3,double,prop> particles_bulk(particles,0);
  vector_dist_subset<3,double,prop> particles_boundary(particles,1);
  auto &bulkIds=particles_bulk.getIds();
  auto &bdrIds=particles_boundary.getIds();

  auto DErr=getV<ERR>(particles);
  auto Af=getV<ANADF>(particles);
  auto Df=getV<DF>(particles);
  auto f=getV<F>(particles);
  auto N=getV<NORMAL>(particles);

  auto verletList = particles.template getVerlet<VL_NON_SYMMETRIC|VL_SKIP_REF_PART>(rCut);
  SurfaceDerivative_x<NORMAL,decltype(verletList)> Sdx{particles,verletList,2,rCut,grid_spacing_surf,rCut/grid_spacing_surf};
  SurfaceDerivative_y<NORMAL,decltype(verletList)> Sdy{particles,verletList,2,rCut,grid_spacing_surf,rCut/grid_spacing_surf};
  SurfaceDerivative_z<NORMAL,decltype(verletList)> Sdz{particles,verletList,2,rCut,grid_spacing_surf,rCut/grid_spacing_surf};

  SurfaceDerivative_xx<NORMAL,decltype(verletList)> Sdxx{particles,verletList,2,rCut,grid_spacing_surf,rCut/grid_spacing_surf};
  SurfaceDerivative_yy<NORMAL,decltype(verletList)> Sdyy{particles,verletList,2,rCut,grid_spacing_surf,rCut/grid_spacing_surf};
  SurfaceDerivative_zz<NORMAL,decltype(verletList)> Sdzz{particles,verletList,2,rCut,grid_spacing_surf,rCut/grid_spacing_surf};

  DErr=-0.5*(Sdx(N[0])+Sdy(N[1])+Sdz(N[2]));
  particles.deleteGhost();
  particles.write_frame(p_output,0,BINARY);


  DCPSE_scheme<equations3d1, decltype(particles)> Solver(particles);
  auto Eig = Sdxx(f)+Sdyy(f)+Sdzz(f)-EIGENVAL*f;
  Solver.impose(Eig, bulkIds, 0);
  Solver.impose(f, bdrIds, 0.1);
  Solver.solve(f);
  particles.deleteGhost();
  particles.write(p_output,BINARY);
  Sdzz.deallocate(particles);
  Sdyy.deallocate(particles);
  Sdxx.deallocate(particles);
  Sdz.deallocate(particles);
  Sdy.deallocate(particles);
  Sdx.deallocate(particles);

  tt.stop();
  if (v_cl.rank() == 0)
    std::cout << "Simulation took: " << tt.getcputime() << " s (CPU) - " << tt.getwct() << " s (wall)\n";
  
  openfpm_finalize();
  return 0;

}