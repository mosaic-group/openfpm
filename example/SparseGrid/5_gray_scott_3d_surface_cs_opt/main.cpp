#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

/*!
 *
 */

constexpr int U = 0;
constexpr int V = 1;
constexpr int phi = 2;
constexpr int normal = 3;
constexpr int tgrad_u = 4;
constexpr int tgrad_v = 5;
constexpr int U_next = 6;
constexpr int V_next = 7;

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

typedef sgrid_dist_soa<3,double,aggregate<double,double,double,double[3],double[3],double[3],double,double> > sgrid_type;

void init(sgrid_type & grid, Box<3,double> & domain)
{
	//! \cond [init sphere channel] \endcond

	auto it = grid.getGridIterator();
	Point<3,double> p1({0.5,0.5,0.5});

	double sx = grid.spacing(0);

	Sphere<3,double> sph1(p1,0.3);
	Sphere<3,double> sph2(p1,0.3 - sx*10);
	Sphere<3,double> sph_zero(p1,0.3 - sx*5);

	while (it.isNext())
	{
		// Get the local grid key
		auto key = it.get_dist();
		auto keyg = it.get();

		Point<3,double> pc;
		Point<3,double> vp;

		for (int i = 0 ; i < 3 ; i++)
        {pc.get(i) = keyg.get(i) * it.getSpacing(i);}

		// Check if the point is in the first sphere
		if (sph1.isInside(pc) == true && sph2.isInside(pc) == false)
		{
			Point<3,double> pn = pc - p1;
			pn /= pn.norm();
			double theta = acos(pn * Point<3,double>({0.0,0.0,1.0}));
			Point<3,double> pn_ = pn;
			pn_[2] = 0.0;
			pn_ /= pn_.norm();
			double aphi = acos(pn_ * Point<3,double>({1.0,0.0,0.0}));

			// Create a perturbation in the solid angle
			if (theta > 0.6 && theta < 0.8 && aphi > 0.0 && aphi < 0.2)
			{
				grid.template insert<U>(key) = 0.5;
				grid.template insert<V>(key) = 0.25;
			}
			else
			{
				grid.template insert<U>(key) = 1.0;
				grid.template insert<V>(key) = 0.0;
			}
			grid.template insert<phi>(key) = sph_zero.distance(pc);
			grid.template insert<normal>(key)[0] = pn[0];
			grid.template insert<normal>(key)[1] = pn[1];
			grid.template insert<normal>(key)[2] = pn[2];

			// Old values U and V
			grid.template insert<U_next>(key) = 0.0;
			grid.template insert<V_next>(key) = 0.0;
		}

		++it;
	}

	//! \cond [init sphere channel] \endcond
}

template<unsigned int U_src,unsigned int V_src,unsigned int U_dst, unsigned int V_dst>
void extend(sgrid_type & grid, size_t (& sz)[3],double (& spacing)[3])
{
	double delta = 1e-10;
	double max = 0.0;

	auto func_extend = [delta,&spacing](auto & grid, auto & ids,
	                                 unsigned char * mask_sum)
									 {
										Vc::double_v phi_c;
										Vc::double_v s;

										Vc::double_v Uext = 0.0;
										Vc::double_v Vext = 0.0;

										Vc::double_v n[3];
										Vc::double_v dir;

										Vc::double_v Uc;
										Vc::double_v Vc;
										Vc::double_v Uc_xm;
										Vc::double_v Vc_xm;
										Vc::double_v Uc_ym;
										Vc::double_v Vc_ym;
										Vc::double_v Uc_zm;
										Vc::double_v Vc_zm;

										Vc::double_v Uc_xp;
										Vc::double_v Vc_xp;
										Vc::double_v Uc_yp;
										Vc::double_v Vc_yp;
										Vc::double_v Uc_zp;
										Vc::double_v Vc_zp;

										load_crs<x,0,phi>(phi_c,grid,ids);
										load_crs_v<x,0,x,normal>(n[x],grid,ids);
										load_crs_v<x,0,y,normal>(n[y],grid,ids);
										load_crs_v<x,0,z,normal>(n[z],grid,ids);

										load_crs<x,0,U_src>(Uc,grid,ids);
										load_crs<x,0,V_src>(Vc,grid,ids);
										load_crs<x,-1,U_src>(Uc_xm,grid,ids);
										load_crs<x,-1,V_src>(Vc_xm,grid,ids);
										load_crs<y,-1,U_src>(Uc_ym,grid,ids);
										load_crs<y,-1,V_src>(Vc_ym,grid,ids);
										load_crs<z,-1,U_src>(Uc_zm,grid,ids);
										load_crs<z,-1,V_src>(Vc_zm,grid,ids);
										load_crs<x,1,U_src>(Uc_xp,grid,ids);
										load_crs<x,1,V_src>(Vc_xp,grid,ids);
										load_crs<y,1,U_src>(Uc_yp,grid,ids);
										load_crs<y,1,V_src>(Vc_yp,grid,ids);
										load_crs<z,1,U_src>(Uc_zp,grid,ids);
										load_crs<z,1,V_src>(Vc_zp,grid,ids);

										s = phi_c / sqrt(phi_c*phi_c + delta*delta);

										dir = s*n[0];
										auto dir_pos = dir > 0;
										auto dir_neg = dir < 0;

										Uext += Vc::iif(dir_pos,dir * (Uc - Uc_xm)/spacing[0],Vc::double_v(0.0));
										Vext += Vc::iif(dir_pos,dir * (Vc - Vc_xm)/spacing[0],Vc::double_v(0.0));
										Uext += Vc::iif(dir_neg,dir * (Uc_xp - Uc)/spacing[0],Vc::double_v(0.0));
										Vext += Vc::iif(dir_neg,dir * (Vc_xp - Vc)/spacing[0],Vc::double_v(0.0));

										dir = s*n[1];
										dir_pos = dir > 0;
										dir_neg = dir < 0;

										Uext += Vc::iif(dir_pos,dir * (Uc - Uc_ym)/spacing[1],Vc::double_v(0.0));
										Vext += Vc::iif(dir_pos,dir * (Vc - Vc_ym)/spacing[1],Vc::double_v(0.0));
										Uext += Vc::iif(dir_neg,dir * (Uc_yp - Uc)/spacing[1],Vc::double_v(0.0));
										Vext += Vc::iif(dir_neg,dir * (Vc_yp - Vc)/spacing[1],Vc::double_v(0.0));

										dir = s*n[2];
										dir_pos = dir > 0;
										dir_neg = dir < 0;

										Uext += Vc::iif(dir_pos,dir * (Uc - Uc_zm)/spacing[2],Vc::double_v(0.0));
										Vext += Vc::iif(dir_pos,dir * (Vc - Vc_zm)/spacing[2],Vc::double_v(0.0));
										Uext += Vc::iif(dir_neg,dir * (Uc_zp - Uc)/spacing[2],Vc::double_v(0.0));
										Vext += Vc::iif(dir_neg,dir * (Vc_zp - Vc)/spacing[2],Vc::double_v(0.0));

										Uext = Uc - 0.0003*Uext;
										Vext = Vc - 0.0003*Vext;

										store_crs<U_dst>(grid,Uext,ids);
										store_crs<V_dst>(grid,Vext,ids);
									 };

		grid.template conv_cross_ids<1,double>({0,0,0},
			{static_cast<long int>(sz[0] - 1),
			 static_cast<long int>(sz[1] - 1),
			 static_cast<long int>(sz[2] - 1)},func_extend);
}

int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	// domain
	Box<3,double> domain({0.0,0.0,0.0},{2.5,2.5,2.5});
	
	// grid size
    size_t sz[3] = {512,512,512};

	// Define periodicity of the grid
	periodicity<3> bc = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
	
	// Ghost in grid unit
	Ghost<3,long int> g(1);
	
	// deltaT
	double deltaT = 0.3;

	// Diffusion constant for specie U
	double du = 1*1e-5;

	// Diffusion constant for specie V
	double dv = 0.5*1e-5;

#ifdef TEST_RUN
        // Number of timesteps
        size_t timeSteps = 200;
#else
	// Number of timesteps
        size_t timeSteps = 100000;
#endif

	// K and F (Physical constant in the equation)
        double K = 0.053;
        double F = 0.014;

	sgrid_type grid(sz,domain,g,bc);

	
	// spacing of the grid on x and y
	double spacing[3] = {grid.spacing(0),grid.spacing(1),grid.spacing(2)};

	init(grid,domain);

	// sync the ghost
	size_t count = 0;
	grid.template ghost_get<U,V>();

	// because we assume that spacing[x] == spacing[y] we use formula 2
	// and we calculate the prefactor of Eq 2
	double uFactor = deltaT * du;
	double vFactor = deltaT * dv;

	auto & v_cl = create_vcluster();

	timer tot_sim;
	tot_sim.start();

	for (size_t i = 0; i < timeSteps ; ++i)
	{
		auto func_grad = [&spacing](auto & grid, auto & ids,
                                 unsigned char * mask_sum){

												Vc::double_v n[3];

												Vc::double_v Uc;
												Vc::double_v xmU;
												Vc::double_v ymU;
												Vc::double_v zmU;

												Vc::double_v Vc;
												Vc::double_v xmV;
												Vc::double_v ymV;
												Vc::double_v zmV;

												Vc::double_v u_out[3];
												Vc::double_v v_out[3];

												load_crs<x,-1,U>(xmU,grid,ids);
												load_crs<y,-1,U>(ymU,grid,ids);
												load_crs<z,-1,U>(zmU,grid,ids);
												load_crs<x,0,U>(Uc,grid,ids);

												load_crs<x,-1,V>(xmV,grid,ids);
												load_crs<y,-1,V>(ymV,grid,ids);
												load_crs<z,-1,V>(zmV,grid,ids);
												load_crs<x,0,V>(Vc,grid,ids);

												load_crs_v<x,0,x,normal>(n[x],grid,ids);
												load_crs_v<x,0,y,normal>(n[y],grid,ids);
												load_crs_v<x,0,z,normal>(n[z],grid,ids);

												u_out[0] = (1.0-n[0]*n[0])*(Uc - xmU)/spacing[0]  + (-n[1]*n[1])*(Uc - ymU)/spacing[1]    + (-n[2]*n[2])*(Uc - zmU)/spacing[2];
												u_out[1] = (-n[0]*n[0])*(Uc - xmU)/spacing[0]    + (1.0-n[1]*n[1])*(Uc - ymU)/spacing[1] + (-n[2]*n[2])*(Uc - zmU)/spacing[2];
												u_out[2] = (-n[0]*n[0])*(Uc - xmU)/spacing[0]    + (-n[1]*n[1])*(Uc - ymU)/spacing[1]    + (1.0-n[2]*n[2])*(Uc - zmU)/spacing[2];

												v_out[0] = (1.0-n[0]*n[0])*(Vc - xmV)/spacing[0] + (-n[1]*n[1])*(Vc - ymV)/spacing[1]    + (-n[2]*n[2])*(Vc - zmV)/spacing[2];
												v_out[1] = (-n[0]*n[0])*(Vc - xmV)/spacing[0]    + (1.0-n[1]*n[1])*(Vc - ymV)/spacing[1] + (-n[2]*n[2])*(Vc - zmV)/spacing[2];
												v_out[2] = (-n[0]*n[0])*(Vc - xmV)/spacing[0]    + (-n[1]*n[1])*(Vc - ymV)/spacing[1]    + (1.0-n[2]*n[2])*(Vc - zmV)/spacing[2];

												Vc::Mask<double> surround;

												for (int i = 0 ; i < Vc::double_v::Size ; i++)
												{surround[i] = (mask_sum[i] == 6);}

												u_out[0] = Vc::iif(surround,u_out[0],Vc::double_v(0.0));
												u_out[1] = Vc::iif(surround,u_out[1],Vc::double_v(0.0));
												u_out[2] = Vc::iif(surround,u_out[2],Vc::double_v(0.0));

												v_out[0] = Vc::iif(surround,v_out[0],Vc::double_v(0.0));
												v_out[1] = Vc::iif(surround,v_out[1],Vc::double_v(0.0));
												v_out[2] = Vc::iif(surround,v_out[2],Vc::double_v(0.0));

												store_crs_v<tgrad_u,x>(grid,u_out[0],ids);
												store_crs_v<tgrad_u,y>(grid,u_out[1],ids);
												store_crs_v<tgrad_u,z>(grid,u_out[2],ids);

												store_crs_v<tgrad_v,x>(grid,v_out[0],ids);
												store_crs_v<tgrad_v,y>(grid,v_out[1],ids);
												store_crs_v<tgrad_v,z>(grid,v_out[2],ids);
											};

		grid.template conv_cross_ids<1,double>({0,0,0},
										  {static_cast<long int>(sz[0] - 1),
										   static_cast<long int>(sz[1] - 1),
										   static_cast<long int>(sz[2] - 1)},func_grad);

		auto func_lap = [&spacing,uFactor,vFactor,deltaT,K,F](auto & grid, auto & ids,
                                 unsigned char * mask_sum){

												Vc::double_v gradU_px;
												Vc::double_v gradU_py;
												Vc::double_v gradU_pz;

												Vc::double_v gradU_x;
												Vc::double_v gradU_y;
												Vc::double_v gradU_z;

												Vc::double_v gradV_px;
												Vc::double_v gradV_py;
												Vc::double_v gradV_pz;

												Vc::double_v gradV_x;
												Vc::double_v gradV_y;
												Vc::double_v gradV_z;

												Vc::double_v lapU;
												Vc::double_v lapV;

												Vc::double_v Uc;
												Vc::double_v Vc;

												Vc::double_v outU;
												Vc::double_v outV;

												load_crs_v<x,1,x,tgrad_u>(gradU_px,grid,ids);
												load_crs_v<y,1,y,tgrad_u>(gradU_py,grid,ids);
												load_crs_v<z,1,z,tgrad_u>(gradU_pz,grid,ids);

												load_crs_v<x,0,x,tgrad_u>(gradU_x,grid,ids);
												load_crs_v<x,0,y,tgrad_u>(gradU_y,grid,ids);
												load_crs_v<x,0,z,tgrad_u>(gradU_z,grid,ids);

												load_crs_v<x,1,x,tgrad_v>(gradV_px,grid,ids);
												load_crs_v<y,1,y,tgrad_v>(gradV_py,grid,ids);
												load_crs_v<z,1,z,tgrad_v>(gradV_pz,grid,ids);

												load_crs_v<x,0,x,tgrad_v>(gradV_x,grid,ids);
												load_crs_v<x,0,y,tgrad_v>(gradV_y,grid,ids);
												load_crs_v<x,0,z,tgrad_v>(gradV_z,grid,ids);

												load_crs<x,0,U>(Uc,grid,ids);
												load_crs<x,0,V>(Vc,grid,ids);

												lapU += (gradU_px - gradU_x) / spacing[0];
												lapV += (gradV_px - gradV_x) / spacing[0];
												lapU += (gradU_py - gradU_y) / spacing[1];
												lapV += (gradV_py - gradV_y) / spacing[1];
												lapU += (gradU_pz - gradU_z) / spacing[2];
												lapV += (gradV_pz - gradV_z) / spacing[2];

												// update based on Eq 2
												outU = Uc + uFactor * lapU +
																			- deltaT * Uc * Vc * Vc +
																			- deltaT * F * (Uc - 1.0);


												// update based on Eq 2
												outV = Vc + vFactor * lapV +
																			deltaT * Uc * Vc * Vc +
																			- deltaT * (F+K) * Vc;

												Vc::Mask<double> surround;

												for (int i = 0 ; i < Vc::double_v::Size ; i++)
												{surround[i] = (mask_sum[i] == 6);}


												outU = Vc::iif(surround,outU,Uc);
												outV = Vc::iif(surround,outV,Vc);

												store_crs<U_next>(grid,outU,ids);
												store_crs<V_next>(grid,outV,ids);
											};

		grid.template conv_cross_ids<1,double>({0,0,0},
										  {static_cast<long int>(sz[0] - 1),
										   static_cast<long int>(sz[1] - 1),
										   static_cast<long int>(sz[2] - 1)},func_lap);

//		New.write_frame("update",i);

		// Extend

		if (i % 5 == 0)
		{
		for (int j = 0 ; j < 2 ; j++)
		{
			if (j % 2 == 0)
			{extend<U_next,V_next,U,V>(grid,sz,spacing);}
			else
			{extend<U,V,U_next,V_next>(grid,sz,spacing);}

			// Here we copy New into the old grid in preparation of the new step
			// It would be better to alternate, but using this we can show the usage
			// of the function copy. To note that copy work only on two grid of the same
			// decomposition. If you want to copy also the decomposition, or force to be
			// exactly the same, use Old = New
			//New.copy_sparse(Old);
		}
		}

/*		auto it = grid.getDomainIterator();

		while (it.isNext())
		{
			// center point
			auto Cp = it.get();

			// update based on Eq 2
			grid.insert<U>(Cp) = grid.get<U_next>(Cp);
			grid.insert<V>(Cp) = grid.get<V_next>(Cp);

			++it;
		}*/

		//! \cond [stencil get and use] \endcond

		// After copy we synchronize again the ghost part U and V
		grid.ghost_get<U,V>();

		// Every 500 time step we output the configuration for
		// visualization
		if (i % 500 == 0)
		{
//			grid.save("output_" + std::to_string(count));
			count++;
		}

		if (v_cl.rank() == 0)
		{std::cout << "STEP: " << i  << "   " << std::endl;}
		if (i % 1000 == 0)
		{
			grid.write_frame("out",i);
		}
	}

	tot_sim.stop();
	std::cout << "Total simulation: " << tot_sim.getwct() << std::endl;

	//! \cond [time stepping] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse Gray Scott in 3D
	 *
	 * ## Finalize ##
	 *
	 * Deinitialize the library
	 *
	 * \snippet Grid/3_gray_scott/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse Gray Scott in 3D
	 *
	 * # Full code # {#code}
	 *
	 * \include Grid/3_gray_scott_3d/main.cpp
	 *
	 */
}
