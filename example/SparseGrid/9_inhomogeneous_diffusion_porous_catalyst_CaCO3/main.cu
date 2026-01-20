//
// Created by jstark on 2023-04-24.
//
/**
 * \file inhomogeneous_diffusion_porous_catalyst_CaCO3/main.cu
 * \page inhomogeneous_diffusion_porous_catalyst_CaCO3 Inhomogeneous diffusion in a CaCO3 particle packed bed
 * [TOC]
 * # Simulate inhomogeneous diffusion in a CaCO\f$_3\f$ particle packed bed #
 *
 * In this example, we simulate inhomogeneous diffusion in the fluid phase of a CaCO\f$_3\f$ particle packed bed. For the complete image-based simulation pipeline and performance of this simulation, we refer to
 * <a href="https://doi.org/10.1016/j.jocs.2023.102118">J. Stark, I. F. Sbalzarini "An open-source pipeline for solving continuous reaction-diffusion models in image-based geometries of porous media", Journal of Computational Science (2023).</a>

 * The geometry of the fluid phase is reconstructed based on \textmuCT images kindly provided by Prof. Jörg Petrasch (Michigan State University, College of Engineering) (see: <a href="https://asmedigitalcollection.asme.org/heattransfer/article/131/7/072701/444766/Tomographic-Characterization-of-a-Semitransparent">S. Haussener et al., “Tomographic characterization of a semitransparent-particle packed bed and determination of its thermal radiative properties”. In: Journal of Heat Transfer 131.7 (2009)</a>). The folder containing the input files for this example (3D dense grid with level-set signed-distance function of the particle packed bed's fluid phase plus textfiles to initialize grid size and dimensions) can be downloaded from our <a href="https://git.mpi-cbg.de/mosaic/software/image-based-modeling/reactiondiffusion_imagebased_porousmedia/-/tree/main/porous_catalysts/data/output_sussman_CaCO3_fluidPhase_531x531x531">git-repository</a>. Please adapt the INPUT_PATH in this code example accordingly.

 * For image based reconstruction and redistancing, see @ref example_sussman_images_3D.
 *
 */

/**
 * \page inhomogeneous_diffusion_porous_catalyst_CaCO3 Inhomogeneous diffusion in a CaCO3 particle packed bed
 *
 * ## Include ## {#e_CaC03_include}
 *
 * We start with inluding header files, setting some paths and indices:
 *
 * \snippet  SparseGrid/inhomogeneous_diffusion_porous_catalyst_CaCO3/main.cu Include
 *
 */
//! \cond [Include] \endcond


#include <iostream>
#include <typeinfo>

#include "CSVReader/CSVReader.hpp"

// Include redistancing files
#include "util/PathsAndFiles.hpp"
#include "level_set/redistancing_Sussman/RedistancingSussman.hpp"
#include "RawReader/InitGridWithPixel.hpp"
#include "include/RemoveLines.hpp" // For removing thin (diagonal or straight) lines

#include "level_set/redistancing_Sussman/HelpFunctionsForGrid.hpp"
#include "level_set/redistancing_Sussman/AnalyticalSDF.hpp"
#include "FiniteDifference/FD_simple.hpp"
#include "Decomposition/Distribution/BoxDistribution.hpp"
#include "timer.hpp"

// Help functions for simulating diffusion in the image-based geometry
#include "include/DiffusionSpace_sparseGrid.hpp"
#include "include/HelpFunctions_diffusion.hpp"	



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// input
const std::string path_to_redistancing_result =
"output_sussman_maxIter6e3_CaCO3_fluidPhase_531x531x531";


const std::string redistancing_filename = "grid_CaCO3_post_redistancing.hdf5";
const std::string path_to_size = path_to_redistancing_result;

// output
const std::string output_name = "output_inhomogDiffusion_CaCO3_fluid_phase";
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Property indices full grid
const size_t PHI_FULL 		 = 0;

typedef aggregate<float> props_full;
const openfpm::vector<std::string> prop_names_full = {"Phi_Sussman_Out"};


// Property indices sparse grid
const size_t PHI_PHASE           = 0;
const size_t CONC_N              = 1;
const size_t CONC_NPLUS1         = 2;
const size_t DIFFUSION_COEFFICIENT  = 3;


typedef aggregate<float, float, float, float> props_sparse;
const openfpm::vector<std::string> prop_names_sparse = {"PHI_PHASE",
                                                        "CONC_N",
                                                        "CONC_NPLUS1",
                                                        "DIFFUSION_COEFFICIENT"};

// Space indices
constexpr size_t x = 0, y = 1, z = 2;

// Parameters for grid size
const size_t dims     = 3;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! \cond [Include] \endcond

/**
 * \page inhomogeneous_diffusion_porous_catalyst_CaCO3 Inhomogeneous diffusion in a CaCO3 particle packed bed
 *
 * ## Initialization and output folder ## {#e_CaC03_init}
 *
 * We start with
 * * Initializing OpenFPM
 * * Setting the output path and creating an output folder
 *
 * \snippet  SparseGrid/inhomogeneous_diffusion_porous_catalyst_CaCO3/main.cu Initialization and output folder
 *
 */
//! \cond [Initialization and output folder] \endcond
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
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Set current working directory, define output paths and create folders where output will be saved
        std::string cwd                     = get_cwd();
        const std::string path_output       = cwd + "/" + output_name + "/";
        create_directory_if_not_exist(path_output);

        if(v_cl.rank()==0) std::cout << "Redistancing result will be loaded from " << path_to_redistancing_result << redistancing_filename << std::endl;

        if(v_cl.rank()==0) std::cout << "Outputs will be saved to " << path_output << std::endl;

	//! \cond [Initialization and output folder] \endcond
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * \page inhomogeneous_diffusion_porous_catalyst_CaCO3 Inhomogeneous diffusion in a CaCO3 particle packed bed
	 *
	 * ## Create a dense grid and load redistancing result ## {#e_CaC03_grid}
	 *
	 * \snippet  SparseGrid/inhomogeneous_diffusion_porous_catalyst_CaCO3/main.cu Grid creation
	 *
	 */
	//! \cond [Grid creation] \endcond
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

        g_dist.load(path_to_redistancing_result + "/" + redistancing_filename); // Load SDF
	//! \cond [Grid creation] \endcond
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * \page inhomogeneous_diffusion_porous_catalyst_CaCO3 Inhomogeneous diffusion in a CaCO3 particle packed bed
	 *
	 * ## Create a sparse grid of fluid phase ## {#e_CaC03_sparse_grid}
	 * 
	 * Obtain sparse grid by placing grid points inside the fluid phase only using the signed distance function
	 *
	 * \snippet  SparseGrid/inhomogeneous_diffusion_porous_catalyst_CaCO3/main.cu Sparse grid creation
	 *
	 */
	//! \cond [Sparse grid creation] \endcond
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Get sparse grid of diffusion domain
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Create sparse grid
        std::cout << "Rank " << v_cl.rank() <<  " starts creating sparse grid." << std::endl;
	typedef sgrid_dist_id_gpu<dims, float, props_sparse, CudaMemory, Dec> sparse_grid_type;
        sparse_grid_type g_sparse(sz, box, ghost);
        g_sparse.setPropNames(prop_names_sparse);

        // Lower and upper bound for SDF indicating fluid phase
        const float d_low = 0.0, d_up = Lx_up;

        // Alternatively: load sparse grid of diffusion domain from HDF5 file
//      g_sparse.load(path_to_redistancing_result + "/" + redistancing_filename);

        // Obtain sparse grid
        get_diffusion_domain_sparse_grid<PHI_FULL, PHI_PHASE>(g_dist, g_sparse, d_low, d_up);
        std::cout << "Rank " << v_cl.rank() <<  " finished creating sparse grid." << std::endl;

	//! \cond [Sparse grid creation] \endcond

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * \page inhomogeneous_diffusion_porous_catalyst_CaCO3 Inhomogeneous diffusion in a CaCO3 particle packed bed
	 *
	 * ## Defining the inhomogeneous diffusion coefficient ## {#e_CaC03_diff_coeff}
	 * 
	 * The inhomogeneous diffusion coefficient is smoothed around the interface by a sigmoidal function that depends on
	 * the signed distance function
	 *
	 * \snippet  SparseGrid/inhomogeneous_diffusion_porous_catalyst_CaCO3/main.cu Diffusion coefficient
	 *
	 */
	//! \cond [Diffusion coefficient] \endcond
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialize properties on sparse grid
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//  Diffusion coefficient
	const float min_porosity = 0.0;
	const float max_porosity = 1.0;
	
	float D_molecular = 1.0;
	float D_min = D_molecular * min_porosity;
	float D_max = D_molecular * max_porosity;
	
	float min_sdf_calcite = get_min_val<PHI_PHASE>(g_sparse);
	
	float x_stretch = 4.0 * g_sparse.spacing(x);
	float x_shift   = min_sdf_calcite * x_stretch;
	float y_min     = D_min;
	float y_max     = D_max - D_min;

	const float c0 = 10.0;
	
	auto dom = g_sparse.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		
		Point<dims, float> coords = g_sparse.getPos(key);
		if(coords.get(x) <= 0.05 * Lx_up) 
		{
			g_sparse.template insertFlush<CONC_N>(key) = c0;
		}
		else g_sparse.template insertFlush<CONC_N>(key) = 0.0;
		
		// Compute diffusion coefficient and store on grid
		g_sparse.template insertFlush<DIFFUSION_COEFFICIENT>(key) =
				get_smooth_sigmoidal(
						g_sparse.template getProp<PHI_PHASE>(key) - min_sdf_calcite,
						x_shift,
						x_stretch,
						y_min,
						y_max);
		++dom;
	}

	//! \cond [Diffusion coefficient] \endcond
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * \page inhomogeneous_diffusion_porous_catalyst_CaCO3 Inhomogeneous diffusion in a CaCO3 particle packed bed
	 *
	 * ## Defining the functor for the inhomogeneous diffusion ## {#e_CaC03_functor}
	 * 
	 * The functor contains the forward time central space discretization of the inhomogeneous diffusion equation.
	 * The stencil size is 1.
	 * It will be passed to the GPU and convolved with the sparse grid.
	 *
	 * \snippet  SparseGrid/inhomogeneous_diffusion_porous_catalyst_CaCO3/main.cu Functor
	 *
	 */
	//! \cond [Functor] \endcond
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Defining the functor that will be passed to the GPU to simulate reaction-diffusion
	// GetCpBlockType<typename SGridGpu, unsigned int prp, unsigned int stencil_size>
	typedef typename GetCpBlockType<decltype(g_sparse),0,1>::type CpBlockType;
	
	const float dx = g_sparse.spacing(x), dy = g_sparse.spacing(y), dz = g_sparse.spacing(z);
	// Stability criterion for FTCS
	const float dt = diffusion_time_step(g_sparse, D_max);
	
	std::cout << "dx = " << dx << "dy = " << dy << "dz = " << dz << std::endl;
	std::cout << "dt = " << dt << std::endl;
	
	
	auto func_inhomogDiffusion = [dx, dy, dz, dt, d_low] __device__ (
			float & u_out,          // concentration output
			float & D_out,          // diffusion coefficient output (dummy)
			float & phi_out,        // sdf of domain output (dummy)
			CpBlockType & u,        // concentration input 
			CpBlockType & D,        // diffusion coefficient input (to get also the half step diffusion coefficients)
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
		
		// Diffusion coefficient
		float D_c  = D(i, j, k);
		float D_px = D(i+1, j, k);
		float D_mx = D(i-1, j, k);
		float D_py = D(i, j+1, k);
		float D_my = D(i, j-1, k);
		float D_pz = D(i, j, k+1);
		float D_mz = D(i, j, k-1);		
		
		// Impose no-flux boundary conditions
		if (phi_px <= d_low + std::numeric_limits<float>::epsilon()) {u_px = u_c; D_px = D_c;}
		if (phi_mx <= d_low + std::numeric_limits<float>::epsilon()) {u_mx = u_c; D_mx = D_c;}
		if (phi_py <= d_low + std::numeric_limits<float>::epsilon()) {u_py = u_c; D_py = D_c;}
		if (phi_my <= d_low + std::numeric_limits<float>::epsilon()) {u_my = u_c; D_my = D_c;}
		if (phi_pz <= d_low + std::numeric_limits<float>::epsilon()) {u_pz = u_c; D_pz = D_c;}
		if (phi_mz <= d_low + std::numeric_limits<float>::epsilon()) {u_mz = u_c; D_mz = D_c;}
		
		// Interpolate diffusion constant of half-step neighboring points
		float D_half_px = (D_c + D_px)/2.0;
		float D_half_mx = (D_c + D_mx)/2.0;
		float D_half_py = (D_c + D_py)/2.0;
		float D_half_my = (D_c + D_my)/2.0;
		float D_half_pz = (D_c + D_pz)/2.0;
		float D_half_mz = (D_c + D_mz)/2.0;
		
		// Compute concentration of next time point
		u_out = u_c + dt *
				(1/(dx*dx) * (D_half_px * (u_px - u_c) - D_half_mx * (u_c - u_mx)) +
				 1/(dy*dy) * (D_half_py * (u_py - u_c) - D_half_my * (u_c - u_my)) +
				 1/(dz*dz) * (D_half_pz * (u_pz - u_c) - D_half_mz * (u_c - u_mz)));

		// dummy outs
		D_out 	= D_c;
		phi_out = phi_c;
	};


	//! \cond [Functor] \endcond
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * \page inhomogeneous_diffusion_porous_catalyst_CaCO3 Inhomogeneous diffusion in a CaCO3 particle packed bed
	 *
	 * ## Run inhomogeneous diffusion ## {#e_CaC03_run}
	 * 
	 *
	 * \snippet  SparseGrid/inhomogeneous_diffusion_porous_catalyst_CaCO3/main.cu Run
	 *
	 */
	//! \cond [Run] \endcond
	// Create text file to which wall clock time for iteration incl. ghost communication will be saved
	std::string path_time_filewct = path_output + "/time_wct_incl_ghostCommunication" + std::to_string(v_cl.rank()) + ".txt";
	std::ofstream out_wct_file(path_time_filewct, std::ios_base::app);
	
        // Create text file to which GPU time will be saved
        std::string path_time_filegpu = path_output + "/time_GPU" + std::to_string(v_cl.rank()) + ".txt";
        std::ofstream out_GPUtime_file(path_time_filegpu, std::ios_base::app);

	// Iterative diffusion
	const size_t iterations = 103;
	const size_t number_write_outs = 10;
	const size_t interval_write = iterations / number_write_outs; // Find interval for writing to reach
	size_t iter = 0;
	float t = 0;

	// Copy from host to GPU for simulation
	g_sparse.template hostToDevice<CONC_N, CONC_NPLUS1, DIFFUSION_COEFFICIENT, PHI_PHASE>();
	g_sparse.template ghost_get<DIFFUSION_COEFFICIENT, PHI_PHASE>(RUN_ON_DEVICE | SKIP_LABELLING);
	t_iterative_diffusion_total.start();
	while(iter <= iterations)
	{
		if (iter % 2 == 0)
		{
			t_iteration_wct.start();
			g_sparse.template ghost_get<CONC_N>(RUN_ON_DEVICE  | SKIP_LABELLING);
			t_GPU.startGPU();
			g_sparse.template conv3_b<
			        CONC_N,
					DIFFUSION_COEFFICIENT,
					PHI_PHASE,
					CONC_NPLUS1,
					DIFFUSION_COEFFICIENT,
					PHI_PHASE, 1>
					({0, 0, 0}, 
					{(long int) sz[x]-1, (long int) sz[y]-1, (long int) sz[z]-1}, 
					func_inhomogDiffusion);
			t_GPU.stopGPU();
			t_iteration_wct.stop();
			
			// Write out time to text-file
			out_wct_file << t_iteration_wct.getwct() << ",";
			out_GPUtime_file << t_GPU.getwctGPU() << ",";
		}
		else
		{
			t_iteration_wct.start();
			g_sparse.template ghost_get<CONC_NPLUS1>(RUN_ON_DEVICE  | SKIP_LABELLING);
			t_GPU.startGPU();
			g_sparse.template conv3_b<
					CONC_NPLUS1,
					DIFFUSION_COEFFICIENT,
					PHI_PHASE,
					CONC_N,
					DIFFUSION_COEFFICIENT,
					PHI_PHASE, 1>
					({0, 0, 0}, 
					{(long int) sz[x]-1, (long int) sz[y]-1, (long int) sz[z]-1}, 
					func_inhomogDiffusion);
                        t_GPU.stopGPU();
                        t_iteration_wct.stop();
			// Write out time to text-file
			out_wct_file << t_iteration_wct.getwct() << ",";
			out_GPUtime_file << t_GPU.getwctGPU() << ",";
		}


		if (iter % interval_write == 0)
		{
			// Copy from GPU to host for writing
			g_sparse.template deviceToHost<CONC_N, CONC_NPLUS1, DIFFUSION_COEFFICIENT, PHI_PHASE>();
			
			// Write g_sparse to vtk
			g_sparse.write_frame(path_output + "g_sparse_diffuse", iter, FORMAT_BINARY);

			// Write g_sparse to hdf5
			g_sparse.save(path_output + "g_sparse_diffuse_" + std::to_string(iter) + ".hdf5");
			
			if (iter % 2 == 0)
			{
				monitor_total_concentration<CONC_N>(g_sparse, t, iter, path_output,"total_conc.csv");
			}
			else
			{
				monitor_total_concentration<CONC_NPLUS1>(g_sparse, t, iter, path_output,"total_conc.csv");
			}
		}


		t += dt;
		iter += 1;
	}
	t_iterative_diffusion_total.stop();
	out_GPUtime_file.close();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	t_total.stop();
	std::cout << "Rank " << v_cl.rank() <<  " total time for the whole programm : " << t_total.getwct() << std::endl;

	//! \cond [Run] \endcond

	/**
	 * \page inhomogeneous_diffusion_porous_catalyst_CaCO3 Inhomogeneous diffusion in a CaCO3 particle packed bed
	 *
	 * ## Terminate ## {#e_CaC03_terminate}
	 *
	 * We end with terminating OpenFPM.
	 *
	 * \snippet  SparseGrid/inhomogeneous_diffusion_porous_catalyst_CaCO3/main.cu Terminate
	 */
	//! \cond [Terminate] \endcond
	openfpm_finalize(); // Finalize openFPM library
	return 0;
}
//! \cond [Terminate] \endcond

/**
 * \page inhomogeneous_diffusion_porous_catalyst_CaCO3 Inhomogeneous diffusion in a CaCO3 particle packed bed
 *
 * ## Full code ## {#e_CaC03_full}
 *
 * @include SparseGrid/inhomogeneous_diffusion_porous_catalyst_CaCO3/main.cu
 */
