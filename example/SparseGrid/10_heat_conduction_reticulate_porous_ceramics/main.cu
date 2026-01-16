//
// Created by jstark on 2023-05-05.
//
/*!
 *
 * \page 10_heat_conduction_RPC Heat conduction in reticulate porous ceramics using sparse grids on GPU
 *
 * [TOC]
 *
 * # Solving heat conduction in the image-based geometry of reticulate porous ceramics # {#e10_heat_conduction_RPC_gpu}
 *
 * In this example, we simulate heat conduction in the solid phase of reticulate porous ceramics with heat dissipation at the surface. For the complete image-based simulation pipeline and performance of this simulation, we refer to
 * <a href="https://doi.org/10.1016/j.jocs.2023.102118">J. Stark, I. F. Sbalzarini "An open-source pipeline for solving continuous reaction-diffusion models in image-based geometries of porous media", Journal of Computational Science (2023).</a>
 * The geometry of the solid phase is reconstructed based on $\upmu$CT images provided by kindly provided by Prof. Jörg Petrasch (Michigan State University, College of Engineering) (see: <a href="https://doi.org/10.1111/j.1551-2916.2008.02308.x">J. Petrasch et al., "Tomography-based multiscale analyses of the 3D geometrical morphology of reticulated porous ceramics", Journal of the American Ceramic Society (2008)</a>
 * and <a href="https://doi.org/10.1115/1.4000226">S. Haussener et al., "Tomography-based heat and mass transfer characterization of reticu- late porous ceramics for high-temperature processing", Journal of Heat Transfer (2010)</a>). The folder containing the input files for this example (3D level-set sparse grid representation of the RPC solid phase plus textfiles to initialize grid size and dimensions) can be downloaded from our <a href="https://git.mpi-cbg.de/mosaic/software/image-based-modeling/reactiondiffusion_imagebased_porousmedia/-/tree/main/porous_catalysts/data/output_sussman_sparse_grid_porousCeramics_1216x1016x941">git-repository</a>. Please adapt the INPUT_PATH in this code example accordingly.
 * 
 *
 * For image based reconstruction and redistancing see @ref example_sussman_images_3D.
 * 
 *
 * First, we initialize and generate sparse grid of solid phase. In this example, an indicator funtion (color field) 
 * would also be sufficient for defining the diffusion domain as we do not scale any parameter by their distance to the
 * surface here, as opposed to the inhomogeneous diffusion in the CaCO\f$_3\f$ fluid phase, 
 * see @ref 9_inhomogeneous_diffusion_porous_catalyst_CaCO3. 
 * We anyways load the SDF as this gives us more flexibility when we want to extend the code towards applications 
 * for which we need the distance to the surface.
 *
 * \snippet SparseGrid/10_heat_conduction_reticulate_porous_ceramics/main.cu initialize
 *
 * Define the functor for heat conduction
 * \snippet SparseGrid/10_heat_conduction_reticulate_porous_ceramics/main.cu functor
 *
 * Iterative heat conduction
 * \snippet SparseGrid/10_heat_conduction_reticulate_porous_ceramics/main.cu iteration
 *
 */

//! \cond [initialize] \endcond
#include <iostream>
#include <typeinfo>

#include "CSVReader/CSVReader.hpp"

// Include redistancing files
#include "util/PathsAndFiles.hpp"
#include "level_set/redistancing_Sussman/RedistancingSussman.hpp"
#include "RawReader/InitGridWithPixel.hpp"
#include "level_set/redistancing_Sussman/HelpFunctionsForGrid.hpp" // For the ghost initialization
#include "include/RemoveLines.hpp" // For removing thin (diagonal or straight) lines

#include "FiniteDifference/FD_simple.hpp"

#include "include/DiffusionSpace_sparseGrid.hpp"
#include "include/HelpFunctions_diffusion.hpp"

#include "Decomposition/Distribution/BoxDistribution.hpp"

const size_t dims     = 3;

// Property indices full grid
const size_t PHI_FULL             = 0;

typedef aggregate<float> props_full;

const openfpm::vector<std::string> prop_names_full = {"Phi_Sussman_Out"};


// Property indices sparse grid
const size_t PHI_PHASE  = 0;
const size_t U_N        = 1;
const size_t U_NPLUS1   = 2;
const size_t K_SINK     = 3;


typedef aggregate<float, float, float, float> props_sparse;
const openfpm::vector<std::string> prop_names_sparse = {"PHI_PHASE",
                                                        "U_N",
                                                        "U_NPLUS1",
                                                        "K_SINK"};


// Space indices
constexpr size_t x = 0, y = 1, z = 2;

// input
const std::string path_to_redistancing_result =
        "output_sussman_sparse_grid_porousCeramics_1216x1016x941/";

const std::string redistancing_filename = "sparseGrid_initial.hdf5";
        
const std::string path_to_size = path_to_redistancing_result;

const std::string output_name = "output_heat_conduction_porous_ceramics_1216x1016x941_sparse_grid";
        

int main(int argc, char* argv[])
{
    // Initialize library.
    openfpm_init(&argc, &argv);
    auto & v_cl = create_vcluster();

    timer t_total;
    timer t_iterative_diffusion_total;
    timer t_iteration_wct;
    timer t_GPU;


    t_total.start();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set current working directory, define output paths and create folders where output will be saved
    std::string cwd                     = get_cwd();
    const std::string path_output       = cwd + "/" + output_name + "/";
    create_directory_if_not_exist(path_output);

    if(v_cl.rank()==0) std::cout << "Redistancing result will be loaded from " << path_to_redistancing_result << redistancing_filename << std::endl;

    if(v_cl.rank()==0) std::cout << "Outputs will be saved to " << path_output << std::endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////
    // Create full grid and load redistancing result
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Read size of sussman redistancing output grid from csv files
    openfpm::vector<size_t> v_sz;
    openfpm::vector<float> v_L;

    size_t m, n;
    read_csv_to_vector(path_to_size + "/N.csv", v_sz, m, n);
    read_csv_to_vector(path_to_size + "/L.csv", v_L, m, n);

    const float Lx_low   = 0.0;
    const float Lx_up    = v_L.get(x);
    const float Ly_low   = 0.0;
    const float Ly_up    = v_L.get(y);
    const float Lz_low   = 0.0;
    const float Lz_up    = v_L.get(z);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    size_t sz[dims] = {v_sz.get(x), v_sz.get(y), v_sz.get(z)};
    Box<dims, float> box({Lx_low, Ly_low, Lz_low}, {Lx_up, Ly_up, Lz_up});
    Ghost<dims, long int> ghost(1);

    // Defining the decomposition and the input grid type
    typedef CartDecomposition<dims,float, CudaMemory, memory_traits_inte, BoxDistribution<dims,float> > Dec;
    typedef grid_dist_id<dims, float, props_full, Dec> grid_in_type;
    grid_in_type g_dist(sz, box, ghost);

    std::cout << "Rank " << v_cl.rank() <<  " starts loading redistancing result..." << std::endl;
    g_dist.load(path_to_redistancing_result + "/" + redistancing_filename); // Load SDF
    std::cout << "Rank " << v_cl.rank() <<  " finished loading redistancing result." << std::endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Get sparse grid of heat conduction domain
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Create sparse grid
    std::cout << "Rank " << v_cl.rank() <<  " starts creating sparse grid." << std::endl;
    

    typedef sgrid_dist_id_gpu<dims, float, props_sparse, CudaMemory, Dec> sparse_grid_type;
    sparse_grid_type g_sparse(sz, box, ghost);
    g_sparse.setPropNames(prop_names_sparse);
    

    const float d_low = 0.0, d_up = Lx_up;
    get_diffusion_domain_sparse_grid<PHI_FULL, PHI_PHASE>(g_dist, g_sparse, d_low, d_up);
    
    // Alternatively: load sparse grid of diffusion domain from HDF5 file
    // g_sparse.load(path_to_redistancing_result + "/" + redistancing_filename);
    
    std::cout << "Rank " << v_cl.rank() <<  " finished creating sparse grid." << std::endl;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Define parameters for heat conduction
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    const float u0          = 1.0; // Initial heat from one side
    const float D           = 0.1; // Heat diffusivity
    const float sink_rate   = 0.1; // Rate of heat loss at solid-air interface
    std::cout << "Diffusion coefficient in s/(mm^2): " << D << std::endl;
    std::cout << "Sink rate in 1/s: " << sink_rate << std::endl;

    // Set initial condition: initial heat source is a sphere of radius R at the domain center
    const float center[dims] = {(Lx_up+Lx_low)/(float)2,
                                (Ly_up+Ly_low)/(float)2,
                                (Lz_up+Lz_low)/(float)2};
    const float R = Lz_up / 4.0;

    auto dom = g_sparse.getDomainIterator();
    while(dom.isNext())
    {
        auto key = dom.get();
        
        // Initial condition: heat comes from sphere at center  
        Point<dims, float> coords = g_sparse.getPos(key);

        if(((coords.get(x) - center[x]) * (coords.get(x) - center[x])
          + (coords.get(y) - center[y]) * (coords.get(y) - center[y])
          + (coords.get(z) - center[z]) * (coords.get(z) - center[z])) <= R * R) 
        {
            g_sparse.template insertFlush<U_N>(key) = u0;
        }
        else g_sparse.template insertFlush<U_N>(key) = 0.0;
        
        ++dom;
    }

    //! \cond [initialize] \endcond

    //! \cond [functor] \endcond
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Defining the functor that will be passed to the GPU to simulate reaction-diffusion
    // GetCpBlockType<typename SGridGpu, unsigned int prp, unsigned int stencil_size>
    typedef typename GetCpBlockType<decltype(g_sparse),0,1>::type CpBlockType;
    
    const float dx = g_sparse.spacing(x), dy = g_sparse.spacing(y), dz = g_sparse.spacing(z);
    // Stability criterion for FTCS : dt = 1/4 * 1/(1/dx^2 + 1/dy^2)
    const float dt = diffusion_time_step(g_sparse, D);
    std::cout << "dx, dy, dz = " << g_sparse.spacing(x) << "," << g_sparse.spacing(y) << "," << g_sparse.spacing(z) << std::endl;
    std::cout << "dt = 1/(4D) * 1/(1/dx^2 + 1/dy^2  + 1/dz^2) = " << dt << std::endl;

    auto func_heatDiffusion_withSink = [dx, dy, dz, dt, d_low, D, sink_rate] __device__ (
            float & u_out,          // concentration output
            float & phi_out,        // sdf of domain output (dummy)
            CpBlockType & u,        // concentration input 
            CpBlockType & phi,      // sdf of domain (for boundary conditions)
            auto & block,
            int offset,
            int i,
            int j,
            int k
    )
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Stencil
        // Concentration
        float u_c  = u(i, j, k);
        float u_px = u(i+1, j, k);
        float u_mx = u(i-1, j, k);
        float u_py = u(i, j+1, k);
        float u_my = u(i, j-1, k);
        float u_pz = u(i, j, k+1);
        float u_mz = u(i, j, k-1);
        
        // Signed distance function
        float phi_c  = phi(i, j, k);
        float phi_px = phi(i+1, j, k);
        float phi_mx = phi(i-1, j, k);
        float phi_py = phi(i, j+1, k);
        float phi_my = phi(i, j-1, k);
        float phi_pz = phi(i, j, k+1);
        float phi_mz = phi(i, j, k-1);

            // Get sink term
            float k_sink = 0; // Sink equal to zero in bulk and equal to sink_rate only at the surfaces
    
        // Impose no-flux boundary conditions and sink term at the surface (i.e., when neighbor is outside)
        if (phi_px <= d_low + std::numeric_limits<float>::epsilon()) {u_px = u_c; k_sink = sink_rate;}
        if (phi_mx <= d_low + std::numeric_limits<float>::epsilon()) {u_mx = u_c; k_sink = sink_rate;}
        if (phi_py <= d_low + std::numeric_limits<float>::epsilon()) {u_py = u_c; k_sink = sink_rate;}
        if (phi_my <= d_low + std::numeric_limits<float>::epsilon()) {u_my = u_c; k_sink = sink_rate;}
        if (phi_pz <= d_low + std::numeric_limits<float>::epsilon()) {u_pz = u_c; k_sink = sink_rate;}
        if (phi_mz <= d_low + std::numeric_limits<float>::epsilon()) {u_mz = u_c; k_sink = sink_rate;}
        
        block.template get<K_SINK>()[offset] = k_sink; // Just for checking in paraview
        
        // Compute concentration of next time point
        u_out = u_c + dt * ( D * ((u_px + u_mx - 2 * u_c)/(dx * dx)
                                + (u_py + u_my - 2 * u_c)/(dy * dy)
                                + (u_pz + u_mz - 2 * u_c)/(dz * dz))
                                - k_sink * u_c);


        // dummy outs
        phi_out = phi_c;
    };
    
    //! \cond [functor] \endcond


    //! \cond [iteration] \endcond
    // Create text file to which wall clock time for iteration incl. ghost communication will be saved
    std::string path_time_filewct = path_output + "/time_wct_incl_ghostCommunication" + std::to_string(v_cl.rank()) + ".txt";
    std::ofstream out_wct_file(path_time_filewct, std::ios_base::app);

    // Create text file to which GPU time will be saved
    std::string path_time_filegpu = path_output + "/time_GPU" + std::to_string(v_cl.rank()) + ".txt";
    std::ofstream out_GPUtime_file(path_time_filegpu, std::ios_base::app);

    
    // Iterative diffusion
    // const size_t iterations = 1e6;
    const size_t iterations = 103;
    const size_t number_write_outs = 10;
    const size_t interval_write = iterations / number_write_outs; // Find interval for writing to reach
    size_t iter = 0;
    float t = 0;

    // Copy from host to GPU for simulation
    g_sparse.template hostToDevice<U_N, U_NPLUS1, K_SINK, PHI_PHASE>();
    g_sparse.template ghost_get<PHI_PHASE, K_SINK>(RUN_ON_DEVICE | SKIP_LABELLING);
    t_iterative_diffusion_total.start();
    while(iter <= iterations)
    {
        if (iter % 2 == 0)
        {
            t_iteration_wct.start();
            g_sparse.template ghost_get<U_N>(RUN_ON_DEVICE | SKIP_LABELLING);
            t_GPU.startGPU();
            g_sparse.template conv2_b<
                    U_N,
                    PHI_PHASE,
                    U_NPLUS1,
                    PHI_PHASE, 1>
                    ({0, 0, 0}, 
                    {(long int) sz[x]-1, (long int) sz[y]-1, (long int) sz[z]-1}, 
                    func_heatDiffusion_withSink);
                    t_GPU.stopGPU();
                    t_iteration_wct.stop();
            
            // Write out time to text-file
            out_wct_file << t_iteration_wct.getwct() << ",";
            out_GPUtime_file << t_GPU.getwctGPU() << ",";

            t_iteration_wct.reset();
        }
        else
        {
            t_iteration_wct.start();
            g_sparse.template ghost_get<U_NPLUS1>(RUN_ON_DEVICE | SKIP_LABELLING);
            t_GPU.startGPU();
            g_sparse.template conv2_b<
                    U_NPLUS1,
                    PHI_PHASE,
                    U_N,
                    PHI_PHASE, 1>
                    ({0, 0, 0}, 
                    {(long int) sz[x]-1, (long int) sz[y]-1, (long int) sz[z]-1}, 
                    func_heatDiffusion_withSink);
            t_GPU.stopGPU();
            t_iteration_wct.stop();
            
            // Write out time to text-file
            out_wct_file << t_iteration_wct.getwct() << ",";
            out_GPUtime_file << t_GPU.getwctGPU() << ",";

            t_iteration_wct.reset();
        }

        
        if (iter % interval_write == 0)
        {
            // Copy from GPU to host for writing
            g_sparse.template deviceToHost<U_N, U_NPLUS1, K_SINK, PHI_PHASE>();
            
            // Write g_sparse to vtk
            g_sparse.write_frame(path_output + "g_sparse_diffuse", iter, FORMAT_BINARY);
            
            if (iter % 2 == 0)
            {
                monitor_total_concentration<U_N>(g_sparse, t, iter, path_output,"total_conc.csv");
            }
            else
            {
                monitor_total_concentration<U_NPLUS1>(g_sparse, t, iter, path_output,"total_conc.csv");
            }
        }

        t += dt;
        iter += 1;
    }
    t_iterative_diffusion_total.stop();
    std::cout << "Rank " << v_cl.rank() <<  " total time for iterative diffusion incl. write outs: " << t_iterative_diffusion_total.getwct() << std::endl;
    std::cout << "Rank " << v_cl.rank() <<  " finished diffusion." << std::endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    t_total.stop();
    std::cout << "Rank " << v_cl.rank() <<  " total time for the whole programm : " << t_total.getwct() << std::endl;


    openfpm_finalize();
    return 0;
}

//! \cond [iteration] \endcond
