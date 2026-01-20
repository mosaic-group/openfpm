/* Benchmark: Redistance a sphere which is only known as indicator function (using OpenFPM + Algoim)
 * Date : 21 June 2021
 * Author : sachin
 */
/**
 *  @file 0_sphere_redistancing/main.cpp
 *  @page Examples_ClosestPoint Redistancing Sphere
 * 
 *  @subpage example_closest_point_sphere
 * 
 */

/**
 * @page example_closest_point_sphere Redistancing a sphere given as indicator function
 * 
 * We are given a grid with a spherical surface implicitly represented using an
 * indicator function, 
 * @f[ f_{ijk} = \begin{cases} 
 *                  1, & \text{inside sphere}\\ 
 *                  -1, & \text{outside sphere} 
 * \end{cases} @f]
 *
 */

#include <cmath>
#include <iostream>
#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "VCluster/VCluster.hpp"
#include "Vector/Vector.hpp"
#include "FiniteDifference/util/common.hpp"
#include "FiniteDifference/eq.hpp"
#include "FiniteDifference/FD_op.hpp"
#include "level_set/closest_point/closest_point.hpp"

constexpr int SIM_DIM = 3;
constexpr int POLY_ORDER = 5;
constexpr int SIM_GRID_SIZE = 64;
constexpr double PI = 3.141592653589793;

// Fields - phi, cp
using GridDist = grid_dist_id<SIM_DIM,double,aggregate<double,double[3]>>;

// Grid size on each dimension
const long int sz[SIM_DIM] = {SIM_GRID_SIZE, SIM_GRID_SIZE, SIM_GRID_SIZE};
const size_t szu[SIM_DIM] = {(size_t) sz[0], (size_t) sz[1], (size_t) sz[2]};

// 2D physical domain
Box<SIM_DIM,double> domain({-1.5,-1.5,-1.5},{1.5,1.5,1.5});

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

// Alias for properties on the grid
constexpr int phi = 0;
constexpr int cp = 1;

double nb_gamma = 0.0;


constexpr int narrow_band_half_width = 8;

typedef struct EllipseParameters{
    double origin[3];
    double radiusA;
    double radiusB;
    double radiusC;
    double eccentricity;
} EllipseParams;

// Generate the initial +1, -1 field for a sphere/ellipsoid
template<const unsigned int phi_field, typename grid_type>
void initializeIndicatorFunc(grid_type &gd, const EllipseParams &params)
{
    // Note: Since we use a Non-periodic boundary, ghost_get does not update ghost layers of sim box.
    // Therefore we are initializing the ghost layer with non-zero values.
    auto it = gd.getDomainGhostIterator();
    while(it.isNext())
    {
        auto key = it.get();
        Point<grid_type::dims, double> coords = gd.getPos(key);
        double posx = coords.get(x);
        double posy = coords.get(y);
        double posz = coords.get(z);
        
        double phi_val = 1.0 - sqrt(((posx - params.origin[0])/params.radiusA)*((posx - params.origin[0])/params.radiusA) + ((posy - params.origin[1])/params.radiusB)*((posy - params.origin[1])/params.radiusB) + ((posz - params.origin[2])/params.radiusC)*((posz - params.origin[2])/params.radiusC));
        gd.template get<phi_field>(key) = (phi_val<0)?-1.0:1.0;
        ++it;
    }
}

// Computes the error in redistancing for a unit sphere
void estimateErrorInReinit(GridDist &gd, EllipseParams &params)
{
    double max_phi_err = -1.0;
    // Update the ls_phi field in ghosts
    gd.template ghost_get<phi>(KEEP_PROPERTIES);

    auto it = gd.getDomainIterator();
    while(it.isNext())
    {
        auto key = it.get();
        Point<GridDist::dims, double> coords = gd.getPos(key);
        double posx = coords.get(0);
        double posy = coords.get(1);
        double posz = coords.get(2);

        double exact_phi_val = 1.0 - sqrt(((posx - params.origin[0])/params.radiusA)*((posx - params.origin[0])/params.radiusA) + ((posy - params.origin[1])/params.radiusB)*((posy - params.origin[1])/params.radiusB) + ((posz - params.origin[2])/params.radiusC)*((posz - params.origin[2])/params.radiusC));
        
        double phi_err = std::abs(exact_phi_val - gd.template get<phi>(key)); // / std::abs(exact_phi);

        if(phi_err > max_phi_err)
            max_phi_err = phi_err;
        ++it;
    }
    
    std::cout<<"Error in redistancing = "<<max_phi_err<<"\n";

    return;
}

int main(int argc, char* argv[])
{
	// Initialize the library
	openfpm_init(&argc, &argv);
    // Create VCluster
    Vcluster<> &v_cl = create_vcluster();

	periodicity<SIM_DIM> grid_bc = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
	// Ghost in grid units
	Ghost <SIM_DIM, long int> grid_ghost(2*narrow_band_half_width);
	GridDist gdist(szu, domain, grid_ghost, grid_bc);

    // Set property names
    openfpm::vector < std::string > prop_names;
    prop_names.add("ls_phi");
    prop_names.add("closest_point");
    gdist.setPropNames(prop_names);

    EllipseParams params;
    params.origin[x] = 0.0;
    params.origin[y] = 0.0;
    params.origin[z] = 0.0;
    params.radiusA = 1.0;
    params.radiusB = 1.0;
    params.radiusC = 1.0;

    initializeIndicatorFunc<phi>(gdist, params);
    nb_gamma = narrow_band_half_width * gdist.spacing(0);
 
    std::cout<<"Grid spacing = "<<gdist.spacing(x)<<"\n";

    // gdist.write_frame("output", 0);

    // Estimates the closest point assuming that the local polynomial approximation
    // of the level set phi field has the correct zero even if it is not the right SDF.
    estimateClosestPoint<phi, cp, POLY_ORDER>(gdist, 3.0);

    // Redistance Levelset - This would try to get the SDF based on CP estimate.
    reinitializeLS<phi, cp>(gdist, 3.0);

    // Compute the errors
    estimateErrorInReinit(gdist, params);

    // gdist.write_frame("output", 1);

	openfpm_finalize();
}

/**
 * @page example_closest_point_sphere Redistancing a sphere given as indicator function
 * 
 * ## Full code ## {#e2d_c_full}
 * 
 * @include example/Numerics/Closest_point/0_sphere_redistancing/main.cpp
 */
