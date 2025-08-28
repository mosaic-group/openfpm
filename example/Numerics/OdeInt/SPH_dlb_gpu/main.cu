/*! \page Numerics SPH Dam break simulation with Dynamic load balacing on Multi-GPU
 *
 *
 * [TOC]
 *
 *
 * # SPH with Dynamic load Balancing on GPU # {#SPH_dlb_gpu}
 *
 *
 * This example show the classical SPH Dam break simulation with load balancing and dynamic load balancing. The main difference with
 * \ref SPH_dlb is that here we use GPUs and 1.2 Millions particles.
 *
 * \htmlonly
 * <a href="#" onclick="hide_show('vector-video-3')" >Simulation video 1</a><br>
 * <div style="display:none" id="vector-video-3">
 * <video id="vid3" width="1200" height="576" controls> <source src="http://openfpm.mpi-cbg.de/web/images/examples/7_SPH_dlb/sph_gpu1.mp4" type="video/mp4"></video>
 * </div>
 * <a href="#" onclick="hide_show('vector-video-4')" >Simulation video 2</a><br>
 * <div style="display:none" id="vector-video-4">
 * <video id="vid4" width="1200" height="576" controls> <source src="http://openfpm.mpi-cbg.de/web/images/examples/7_SPH_dlb/sph_gpu2.mp4" type="video/mp4"></video>
 * </div>
 * <a href="#" onclick="hide_show('vector-video-15')" >Simulation video 3</a><br>
 * <div style="display:none" id="vector-video-15">
 * <video id="vid15" width="1200" height="576" controls> <source src="http://openfpm.mpi-cbg.de/web/images/examples/7_SPH_dlb/sph_gpu3.mp4" type="video/mp4"></video>
 * </div>
 * \endhtmlonly
 *
 * This example use all the features explained in example \ref e3_md_gpu. Additionally this example show how to remove particles
 * on GPU using a bulk remove function on GPU
 *
 * ## Bulk remove
 *
 * On SPH we have the necessity to remove particles that go out of bound. OpenFPM provide the function \b remove_marked \b .
 *
 * \snippet example/Numerics/Odeint/SPH_dlb_gpu/main.cu remove_marked_part
 *
 * where vectorDist is the vector_dist_gpu red is the property that mark which particle must be removed. We mark the particle to be removed in the function kernel
 * We check if the particle go out of the region of interest or their density go critically far from the rest density
 *
 * \snippet example/Numerics/Odeint/SPH_dlb_gpu/main.cu mark_to_remove_kernel
 *
 * ## Macro CUDA_LAUNCH
 *
 * When we want to launch a kernel "my_kernel" on CUDA we in general use the Nvidia CUDA syntax
 *
 * my_kernel<<<wthr,thr>>>(arguments ... )
 *
 * Where wthr is the number of workgroups and thr is the number of threads in a workgroup and arguments... are the arguments to pass to the kernel.
 * Equivalently we can launch a kernel with the macro CUDA_LAUNCH_DIM3(my_kernel,wthr,thr,arguments...) or CUDA_LAUNCH(my_kernel,ite,arguments) where
 * ite has been taken using getDomainIteratorGPU. There are several advantage on using CUDA_LAUNCH. The first advantage in using the macro is enabling SE_CLASS1
 * all kernel launch become synchronous and an error check is performed before continue to the next kernel making debugging easier. Another feature is the possibility
 * to run CUDA code on CPU without a GPU. compiling with "CUDA_ON_CPU=1 make" (Note openfpm must be compiled with GPU support (-g)  or with CUDA_ON_CPU support
 * (-c "... --enable_cuda_on_cpu"). You can compile this example on CPU. You do not have to change a single line of code for this example. (Check the video to see this
 * feature in action). All the openfpm GPU example and CUDA example can run on CPU if they use CUDA_LAUNCH as macro. We are planning to support
 * AMD GPUs as well using this system.
 *
 * \include example/Numerics/Odeint/SPH_dlb_gpu/main.cu
 *
 */

#ifdef __NVCC__

#include <math.h>

#include "Vector/vector_dist.hpp"
#include "Draw/DrawParticles.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"

typedef float real_number;

// A constant to indicate boundary particles
#define BOUNDARY 0

// A constant to indicate fluid particles
#define FLUID 1

// initial spacing between particles dp in the formulas
const real_number dp = 0.0085;
// Maximum height of the fluid water
// is going to be calculated and filled later on
real_number h_swl = 0.0;

// c_s in the formulas (constant used to calculate the sound speed)
const real_number coeff_sound = 20.0;

// gamma in the formulas
const real_number gamma_ = 7.0;

// sqrt(3.0*dp*dp) support of the kernel
const real_number H = 0.0147224318643;

// Eta in the formulas
const real_number Eta2 = 0.01 * H*H;

// alpha in the formula
const real_number visco = 0.1;

// cbar in the formula (calculated later)
real_number cbar = 0.0;

// Mass of the fluid particles
const real_number MassFluid = 0.000614125;

// Mass of the boundary particles
const real_number MassBound = 0.000614125;

// End simulation time
#ifdef TEST_RUN
const real_number t_end = 0.001;
#else
const real_number t_end = 1.5;
#endif

// Gravity acceleration
const real_number gravity = 9.81;

// Reference densitu 1000Kg/m^3
const real_number RhoZero = 1000.0;

// Filled later require h_swl, it is b in the formulas
real_number B = 0.0;

// Constant used to define time integration
const real_number CFLnumber = 0.2;

// Minimum T
const real_number DtMin = 0.00001;

// Minimum Rho allowed
const real_number RhoMin = 700.0;

// Maximum Rho allowed
const real_number RhoMax = 1300.0;

// Filled in initialization
real_number max_fluid_height = 0.0;

// Properties

// FLUID or BOUNDARY
const size_t TYPE = 0;

// Density
const int RHO = 1;

// Density at step n-1
const int RHO_PREV = 2;

// Pressure
const int PRESSURE = 3;

// Delta rho calculated in the force calculation
const int DRHO = 4;

// calculated force
const int FORCE = 5;

// velocity
const int VELOCITY = 6;

// velocity at previous step
const int VELOCITY_PREV = 7;

// temporal variable to store velocity
const int VELOCITY_TMP = 11;

// temporal variable to store density
const int RHO_TMP = 10;

const int RED = 8;

const int RED2 = 9;

// Type of the vector containing particles
typedef vector_dist_gpu<3,real_number,aggregate<size_t,real_number,  real_number,    real_number,     real_number,     VectorS<3, real_number>, VectorS<3, real_number>, VectorS<3, real_number>, real_number, real_number, real_number, VectorS<3, real_number>>> particles;
//                                              |          |             |               |                |                      |                         |                        |                  |           |			|				|
//                                              |          |             |               |                |                      |                         |                        |                  |           |			|				|
//                                             type      density       density        Pressure          delta                  force                    velocity                 velocity           reduction    another	 temp density temp velocity
//                                                                     at n-1                           density                                                                  at n - 1           buffer   reduction buffer


struct ModelCustom
{
	template<typename Decomposition, typename vector>
	inline void addComputation(
		Decomposition & dec,
		vector & vectorDist,
		size_t v,
		size_t p)
	{
		if (vectorDist.template getProp<TYPE>(p) == FLUID)
			dec.addComputationCost(v,4);
		else
			dec.addComputationCost(v,3);
	}

	template<typename Decomposition> inline void applyModel(Decomposition & dec, size_t v)
	{
		dec.setSubSubDomainComputationCost(v, dec.getSubSubDomainComputationCost(v) * dec.getSubSubDomainComputationCost(v));
	}

	real_number distributionTol()
	{
		return 1.01;
	}
};

template<typename vd_type>
__global__ void EqState_gpu(vd_type vectorDist, real_number B)
{
	auto a = GET_PARTICLE(vectorDist);

	real_number rho_a = vectorDist.template getProp<RHO>(a);
	real_number rho_frac = rho_a / RhoZero;

	vectorDist.template getProp<PRESSURE>(a) = B*( rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac - 1.0);
}

inline void EqState(particles & vectorDist)
{
	auto it = vectorDist.getDomainIteratorGPU();

	CUDA_LAUNCH(EqState_gpu,it,vectorDist.toKernel(),B);
}


const real_number a2 = 1.0/M_PI/H/H/H;

inline __device__ __host__ real_number Wab(real_number r)
{
	r /= H;

	if (r < 1.0)
		return (1.0 - 3.0/2.0*r*r + 3.0/4.0*r*r*r)*a2;
	else if (r < 2.0)
		return (1.0/4.0*(2.0 - r)*(2.0 - r)*(2.0 - r))*a2;
	else
		return 0.0;
}


const real_number c1 = -3.0/M_PI/H/H/H/H;
const real_number d1 = 9.0/4.0/M_PI/H/H/H/H;
const real_number c2 = -3.0/4.0/M_PI/H/H/H/H;
const real_number a2_4 = 0.25*a2;
// Filled later
real_number W_dap = 0.0;

inline __device__ __host__ void DWab(Point<3,real_number> & dx, Point<3,real_number> & DW, real_number r)
{
	const real_number qq=r/H;

    real_number qq2 = qq * qq;
    real_number fac1 = (c1*qq + d1*qq2)/r;
    real_number b1 = (qq < 1.0f)?1.0f:0.0f;

    real_number wqq = (2.0f - qq);
    real_number fac2 = c2 * wqq * wqq / r;
    real_number b2 = (qq >= 1.0f && qq < 2.0f)?1.0f:0.0f;

    real_number factor = (b1*fac1 + b2*fac2);

    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
    DW.get(2) = factor * dx.get(2);
}

// Tensile correction
inline __device__ __host__  real_number Tensile(real_number r, real_number rhoa, real_number rhob, real_number prs1, real_number prs2, real_number W_dap)
{
	const real_number qq=r/H;
	//-Cubic Spline kernel
	real_number wab;
	if(r>H)
	{
		real_number wqq1=2.0f-qq;
		real_number wqq2=wqq1*wqq1;

		wab=a2_4*(wqq2*wqq1);
	}
	else
	{
	    real_number wqq2=qq*qq;
	    real_number wqq3=wqq2*qq;

	    wab=a2*(1.0f-1.5f*wqq2+0.75f*wqq3);
	}

	//-Tensile correction.
	real_number fab=wab*W_dap;
	fab*=fab; fab*=fab; //fab=fab^4
	const real_number tensilp1=(prs1/(rhoa*rhoa))*(prs1>0.0f? 0.01f: -0.2f);
	const real_number tensilp2=(prs2/(rhob*rhob))*(prs2>0.0f? 0.01f: -0.2f);

	return (fab*(tensilp1+tensilp2));
}


inline __device__ __host__ real_number Pi(const Point<3,real_number> & dr, real_number rr2, Point<3,real_number> & dv, real_number rhoa, real_number rhob, real_number massb, real_number cbar, real_number & visc)
{
	const real_number dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
	const real_number dot_rr2 = dot/(rr2+Eta2);
	visc=(dot_rr2 < visc)?visc:dot_rr2;

	if(dot < 0)
	{
		const float amubar=H*dot_rr2;
		const float robar=(rhoa+rhob)*0.5f;
		const float pi_visc=(-visco*cbar*amubar/robar);

		return pi_visc;
    }
	else
		return 0.0f;
}

template<typename particles_type, typename CellList_type>
__global__ void calc_forces_gpu(particles_type vectorDist, CellList_type cellList, real_number W_dap, real_number cbar)
{
	auto a = GET_PARTICLE(vectorDist);

	real_number max_visc = 0.0f;

	// Get the position xp of the particle
	Point<3,real_number> xa = vectorDist.getPos(a);

	// Type of the particle
	unsigned int typea = vectorDist.template getProp<TYPE>(a);

	// Get the density of the of the particle a
	real_number rhoa = vectorDist.template getProp<RHO>(a);

	// Get the pressure of the particle a
	real_number Pa = vectorDist.template getProp<PRESSURE>(a);

	// Get the Velocity of the particle a
	Point<3,real_number> va = vectorDist.template getProp<VELOCITY>(a);

	// Reset the force counter (- gravity on zeta direction)
	Point<3,real_number> force_;
	force_.get(0) = 0.0f;
	force_.get(1) = 0.0f;
	force_.get(2) = -gravity;
	real_number drho_ = 0.0f;

	// Get an iterator over the neighborhood particles of p
	auto Np = cellList.getNNIteratorBox(cellList.getCell(xa));

	// For each neighborhood particle
	while (Np.isNext() == true)
	{
		// ... q
		auto b = Np.get();

		// Get the position xp of the particle
		Point<3,real_number> xb = vectorDist.getPos(b);

		if (a == b)	{++Np; continue;};

		unsigned int typeb = vectorDist.template getProp<TYPE>(b);

		real_number massb = (typeb == FLUID)?MassFluid:MassBound;
		Point<3,real_number> vb = vectorDist.template getProp<VELOCITY>(b);
		real_number Pb = vectorDist.template getProp<PRESSURE>(b);
		real_number rhob = vectorDist.template getProp<RHO>(b);

		// Get the distance between p and q
		Point<3,real_number> dr = xa - xb;
		// take the norm of this vector
		real_number r2 = norm2(dr);

		// if they interact
		if (r2 < 4.0*H*H && r2 >= 1e-16)
		{
			real_number r = sqrt(r2);

			Point<3,real_number> v_rel = va - vb;

			Point<3,real_number> DW;
			DWab(dr,DW,r);

			real_number factor = - massb*((Pa + Pb) / (rhoa * rhob) + Tensile(r,rhoa,rhob,Pa,Pb,W_dap) + Pi(dr,r2,v_rel,rhoa,rhob,massb,cbar,max_visc));

			// Bound - Bound does not produce any change
			// factor = (typea == BOUNDARY && typeb == BOUNDARY)?0.0f:factor;
			factor = (typea != FLUID)?0.0f:factor;

			force_.get(0) += factor * DW.get(0);
			force_.get(1) += factor * DW.get(1);
			force_.get(2) += factor * DW.get(2);

			real_number scal = massb*(v_rel.get(0)*DW.get(0)+v_rel.get(1)*DW.get(1)+v_rel.get(2)*DW.get(2));
			scal = (typea == BOUNDARY && typeb == BOUNDARY)?0.0f:scal;

			drho_ += scal;
		}

		++Np;
	}

	vectorDist.template getProp<RED>(a) = max_visc;

	vectorDist.template getProp<FORCE>(a)[0] = force_.get(0);
	vectorDist.template getProp<FORCE>(a)[1] = force_.get(1);
	vectorDist.template getProp<FORCE>(a)[2] = force_.get(2);
	vectorDist.template getProp<DRHO>(a) = drho_;
}

template<typename CellList> inline void calc_forces(particles & vectorDist, CellList & cellList, real_number & max_visc, size_t cnt)
{
	auto part = vectorDist.getDomainIteratorGPU(32);

	// Update the cell-list
	vectorDist.updateCellListGPU(cellList);

	CUDA_LAUNCH(calc_forces_gpu,part,vectorDist.toKernel(),cellList.toKernel(),W_dap,cbar);

	max_visc = reduce_local<RED,_max_>(vectorDist);
}

template<typename vector_type>
__global__ void max_acceleration_and_velocity_gpu(vector_type vectorDist)
{
	auto a = GET_PARTICLE(vectorDist);

	Point<3,real_number> acc(vectorDist.template getProp<FORCE>(a));
	vectorDist.template getProp<RED>(a) = norm(acc);

	Point<3,real_number> vel(vectorDist.template getProp<VELOCITY>(a));
	vectorDist.template getProp<RED2>(a) = norm(vel);
}

template<typename vector_type>
__global__ void checkGPU(vector_type vector)
{
	int i = threadIdx.x;
	printf("Check GPU %d %p %f\n", i, &vector.get<0>(i), vector.get<0>(i));
}



void max_acceleration_and_velocity(particles & vectorDist, real_number & max_acc, real_number & max_vel)
{
	// Calculate the maximum acceleration
	auto part = vectorDist.getDomainIteratorGPU();

	CUDA_LAUNCH(max_acceleration_and_velocity_gpu,part,vectorDist.toKernel());

	max_acc = reduce_local<RED,_max_>(vectorDist);
	max_vel = reduce_local<RED2,_max_>(vectorDist);

	Vcluster<> & v_cl = create_vcluster();
	v_cl.max(max_acc);
	v_cl.max(max_vel);
	v_cl.execute();
}


real_number calc_deltaT(particles & vectorDist, real_number ViscDtMax)
{
	real_number Maxacc = 0.0;
	real_number Maxvel = 0.0;
	max_acceleration_and_velocity(vectorDist,Maxacc,Maxvel);

	//-dt1 depends on force per unit mass.
	const real_number dt_f = (Maxacc)?sqrt(H/Maxacc):std::numeric_limits<float>::max();

	//-dt2 combines the Courant and the viscous time-step controls.
	const real_number dt_cv = H/(std::max(cbar,Maxvel*10.f) + H*ViscDtMax);

	//-dt new value of time step.
	real_number dt=real_number(CFLnumber)*std::min(dt_f,dt_cv);
	if(dt<real_number(DtMin))
	{dt=real_number(DtMin);}

	return dt;
}

template<typename vector_dist_type>
__global__ void checkPosPrpLimits(vector_dist_type vectorDist)
{
	auto p = GET_PARTICLE(vectorDist);

	// if the particle type is boundary
	if (vectorDist.template getProp<TYPE>(p) == BOUNDARY)
	{
		real_number rho = vectorDist.template getProp<RHO>(p);
		if (rho < RhoZero)
			vectorDist.template getProp<RHO>(p) = RhoZero;

		return;
	}

	// Check if the particle go out of range in space and in density, if they do mark them to remove it later
	if (vectorDist.getPos(p)[0] <  0.0 || vectorDist.getPos(p)[1] < 0.0 || vectorDist.getPos(p)[2] < 0.0 ||
		vectorDist.getPos(p)[0] >  1.61 || vectorDist.getPos(p)[1] > 0.68 || vectorDist.getPos(p)[2] > 0.50 ||
		vectorDist.template getProp<RHO>(p) < RhoMin || vectorDist.template getProp<RHO>(p) > RhoMax)
		{vectorDist.template getProp<RED>(p) = 1;}
}

size_t cnt = 0;

void verlet_int(particles & vectorDist, real_number dt)
{
	// particle iterator
	auto part = vectorDist.getDomainIteratorGPU();

	real_number dt205 = dt*dt*0.5;
	real_number dt2 = dt*2.0;

	auto posExpression = getV<POS_PROP, comp_dev>(vectorDist);
	auto posExpression2 = getV<POS_PROP>(vectorDist);
	auto forceExpression = getV<FORCE, comp_dev>(vectorDist);
	auto drhoExpression = getV<DRHO, comp_dev>(vectorDist);
	auto typeExpression = getV<TYPE, comp_dev>(vectorDist);
	auto velocityExpression = getV<VELOCITY, comp_dev>(vectorDist);

	auto rho_tmpExpression = getV<RHO_TMP, comp_dev>(vectorDist);
	auto rhoExpression = getV<RHO, comp_dev>(vectorDist);
	auto rho_prevExpression = getV<RHO_PREV, comp_dev>(vectorDist);

	auto velocity_prevExpression = getV<VELOCITY_PREV, comp_dev>(vectorDist);
	auto velocity_tmpExpression = getV<VELOCITY_TMP, comp_dev>(vectorDist);
	auto redExpression = getV<RED, comp_dev>(vectorDist);

	rho_tmpExpression = rhoExpression;
	rhoExpression = rho_prevExpression + dt2*drhoExpression;
	rho_prevExpression = rho_tmpExpression;

	posExpression = posExpression + velocityExpression*dt + forceExpression*dt205 * typeExpression;

	velocity_tmpExpression = velocityExpression;
	velocityExpression = velocity_prevExpression + forceExpression*dt2 * typeExpression;
	velocity_prevExpression = velocity_tmpExpression;

	redExpression = 0;

	CUDA_LAUNCH(checkPosPrpLimits,part,vectorDist.toKernel());

	// remove the particles marked
	remove_marked<RED>(vectorDist);

	// increment the iteration counter
	cnt++;
}


void euler_int(particles & vectorDist, real_number dt)
{

	// particle iterator
	auto part = vectorDist.getDomainIteratorGPU();

	real_number dt205 = dt*dt*0.5;

	auto posExpression = getV<POS_PROP, comp_dev>(vectorDist);
	auto forceExpression = getV<FORCE, comp_dev>(vectorDist);
	auto drhoExpression = getV<DRHO, comp_dev>(vectorDist);
	auto typeExpression = getV<TYPE, comp_dev>(vectorDist);

	auto rhoExpression = getV<RHO, comp_dev>(vectorDist);
	auto rho_prevExpression = getV<RHO_PREV, comp_dev>(vectorDist);

	auto velocityExpression = getV<VELOCITY, comp_dev>(vectorDist);
	auto velocity_prevExpression = getV<VELOCITY_PREV, comp_dev>(vectorDist);

	auto redExpression = getV<RED, comp_dev>(vectorDist);

	rho_prevExpression = rhoExpression;
	rhoExpression = rhoExpression + dt*drhoExpression;

	posExpression = posExpression + velocityExpression*dt + forceExpression*dt205 * typeExpression;

	velocity_prevExpression = velocityExpression;
	velocityExpression = velocityExpression + forceExpression*dt * typeExpression;

	redExpression = 0;

	CUDA_LAUNCH(checkPosPrpLimits,part,vectorDist.toKernel());

	// remove the particles
	remove_marked<RED>(vectorDist);

	cnt++;
}

template<typename vector_type, typename CellList_type>
__global__ void sensor_pressure_gpu(vector_type vectorDist, CellList_type cellList, Point<3,real_number> probe, real_number * press_tmp)
{
	real_number tot_ker = 0.0;

	// Get the position of the probe i
	Point<3,real_number> xp = probe;

	// get the iterator over the neighbohood particles of the probes position
	auto itg = cellList.getNNIteratorBox(cellList.getCell(xp));
	while (itg.isNext())
	{
		auto q = itg.get();

		// Only the fluid particles are importants
		if (vectorDist.template getProp<TYPE>(q) != FLUID)
		{
			++itg;
			continue;
		}

		// Get the position of the neighborhood particle q
		Point<3,real_number> xq = vectorDist.getPos(q);

		// Calculate the contribution of the particle to the pressure
		// of the probe
		real_number r = sqrt(norm2(xp - xq));

		real_number ker = Wab(r) * (MassFluid / RhoZero);

		// Also keep track of the calculation of the summed
		// kernel
		tot_ker += ker;

		// Add the total pressure contribution
		*press_tmp += vectorDist.template getProp<PRESSURE>(q) * ker;

		// next neighborhood particle
		++itg;
	}

	// We calculate the pressure normalizing the
	// sum over all kernels
	if (tot_ker == 0.0)
		{*press_tmp = 0.0;}
	else
		{*press_tmp = 1.0 / tot_ker * *press_tmp;}
}

template<typename Vector, typename CellList>
inline void sensor_pressure(Vector & vectorDist,
	CellList & cellList,
	openfpm::vector<openfpm::vector<real_number>> & press_t,
	openfpm::vector<Point<3,real_number>> & probes)
{
    Vcluster<> & v_cl = create_vcluster();

    press_t.add();

    for (size_t i = 0 ; i < probes.size() ; i++)
    {
    	// A float variable to calculate the pressure of the problem
    	CudaMemory press_tmp_(sizeof(real_number));
    	real_number press_tmp;

        // if the probe is inside the processor domain
		if (vectorDist.getDecomposition().isLocal(probes.get(i)) == true)
		{
			vectorDist.updateCellListGPU(cellList);

			Point<3,real_number> probe = probes.get(i);
			CUDA_LAUNCH_DIM3(sensor_pressure_gpu,1,1,vectorDist.toKernel(),cellList.toKernel(),probe,(real_number *)press_tmp_.toKernel());

			// move calculated pressure on
			press_tmp_.deviceToHost();
			press_tmp = *(real_number *)press_tmp_.getPointer();
		}

		// This is not necessary in principle, but if you
		// want to make all processor aware of the history of the calculated
		// pressure we have to execute this
		v_cl.sum(press_tmp);
		v_cl.execute();

		// We add the calculated pressure into the history
		press_t.last().add(press_tmp);
	}
}

int main(int argc, char* argv[])
{
    // initialize the library
	openfpm_init(&argc,&argv);

	// It contain for each time-step the value detected by the probes
	openfpm::vector<openfpm::vector<real_number>> press_t;
	openfpm::vector<Point<3,real_number>> probes;

	probes.add({0.8779f,0.3f,0.02f});
	probes.add({0.754f,0.31f,0.02f});

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<3,real_number> domain({-0.05f,-0.05f,-0.05f},{1.7010f,0.7065f,0.511f});
	size_t sz[3] = {207,90,66};

	// Fill W_dap
	W_dap = 1.0/Wab(H/1.5);

	// Here we define the boundary conditions of our problem
    size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<3,real_number> g(2*H);

	particles vectorDist(0,domain,bc,g,DEC_GRAN(512));

	//! \cond [draw fluid] \endcond

	// You can ignore all these dp/2.0 is a trick to reach the same initialization
	// of Dual-SPH that use a different criteria to draw particles
	Box<3,real_number> fluid_box({dp/2.0f,dp/2.0f,dp/2.0f},{0.4f+dp/2.0f,0.67f-dp/2.0f,0.3f+dp/2.0f});

	// return an iterator to the fluid particles to add to vectorDist
	auto fluid_it = DrawParticles::DrawBox(vectorDist,sz,domain,fluid_box);

	// here we fill some of the constants needed by the simulation
	max_fluid_height = fluid_it.getBoxMargins().getHigh(2);
	h_swl = fluid_it.getBoxMargins().getHigh(2) - fluid_it.getBoxMargins().getLow(2);
	B = (coeff_sound)*(coeff_sound)*gravity*h_swl*RhoZero / gamma_;
	cbar = coeff_sound * sqrt(gravity * h_swl);

	// for each particle inside the fluid box ...
	while (fluid_it.isNext())
	{
		// ... add a particle ...
		vectorDist.add();

		// ... and set it position ...
		vectorDist.getLastPos()[0] = fluid_it.get().get(0);
		vectorDist.getLastPos()[1] = fluid_it.get().get(1);
		vectorDist.getLastPos()[2] = fluid_it.get().get(2);

		// and its type.
		vectorDist.template getLastProp<TYPE>() = FLUID;

		// We also initialize the density of the particle and the hydro-static pressure given by
		//
		// RhoZero*g*h = P
		//
		// rho_p = (P/B + 1)^(1/Gamma) * RhoZero
		//

		vectorDist.template getLastProp<PRESSURE>() = RhoZero * gravity *  (max_fluid_height - fluid_it.get().get(2));

		vectorDist.template getLastProp<RHO>() = pow(vectorDist.template getLastProp<PRESSURE>() / B + 1, 1.0/gamma_) * RhoZero;
		vectorDist.template getLastProp<RHO_PREV>() = vectorDist.template getLastProp<RHO>();
		vectorDist.template getLastProp<VELOCITY>()[0] = 0.0;
		vectorDist.template getLastProp<VELOCITY>()[1] = 0.0;
		vectorDist.template getLastProp<VELOCITY>()[2] = 0.0;

		vectorDist.template getLastProp<VELOCITY_PREV>()[0] = 0.0;
		vectorDist.template getLastProp<VELOCITY_PREV>()[1] = 0.0;
		vectorDist.template getLastProp<VELOCITY_PREV>()[2] = 0.0;

		// next fluid particle
		++fluid_it;
	}

	// Recipient
	Box<3,real_number> recipient1({0.0f,0.0f,0.0f},{1.6f+dp/2.0f,0.67f+dp/2.0f,0.4f+dp/2.0f});
	Box<3,real_number> recipient2({dp,dp,dp},{1.6f-dp/2.0f,0.67f-dp/2.0f,0.4f+dp/2.0f});

	Box<3,real_number> obstacle1({0.9f,0.24f-dp/2.0f,0.0f},{1.02f+dp/2.0f,0.36f,0.45f+dp/2.0f});
	Box<3,real_number> obstacle2({0.9f+dp,0.24f+dp/2.0f,0.0f},{1.02f-dp/2.0f,0.36f-dp,0.45f-dp/2.0f});
	Box<3,real_number> obstacle3({0.9f+dp,0.24f,0.0f},{1.02f,0.36f,0.45f});

	openfpm::vector<Box<3,real_number>> holes;
	holes.add(recipient2);
	holes.add(obstacle1);
	auto bound_box = DrawParticles::DrawSkin(vectorDist,sz,domain,holes,recipient1);

	while (bound_box.isNext())
	{
		vectorDist.add();

		vectorDist.getLastPos()[0] = bound_box.get().get(0);
		vectorDist.getLastPos()[1] = bound_box.get().get(1);
		vectorDist.getLastPos()[2] = bound_box.get().get(2);

		vectorDist.template getLastProp<TYPE>() = BOUNDARY;
		vectorDist.template getLastProp<RHO>() = RhoZero;
		vectorDist.template getLastProp<RHO_PREV>() = RhoZero;
		vectorDist.template getLastProp<VELOCITY>()[0] = 0.0;
		vectorDist.template getLastProp<VELOCITY>()[1] = 0.0;
		vectorDist.template getLastProp<VELOCITY>()[2] = 0.0;

		vectorDist.template getLastProp<VELOCITY_PREV>()[0] = 0.0;
		vectorDist.template getLastProp<VELOCITY_PREV>()[1] = 0.0;
		vectorDist.template getLastProp<VELOCITY_PREV>()[2] = 0.0;

		++bound_box;
	}

	auto obstacle_box = DrawParticles::DrawSkin(vectorDist,sz,domain,obstacle2,obstacle1);

	while (obstacle_box.isNext())
	{
		vectorDist.add();

		vectorDist.getLastPos()[0] = obstacle_box.get().get(0);
		vectorDist.getLastPos()[1] = obstacle_box.get().get(1);
		vectorDist.getLastPos()[2] = obstacle_box.get().get(2);

		vectorDist.template getLastProp<TYPE>() = BOUNDARY;
		vectorDist.template getLastProp<RHO>() = RhoZero;
		vectorDist.template getLastProp<RHO_PREV>() = RhoZero;
		vectorDist.template getLastProp<VELOCITY>()[0] = 0.0;
		vectorDist.template getLastProp<VELOCITY>()[1] = 0.0;
		vectorDist.template getLastProp<VELOCITY>()[2] = 0.0;

		vectorDist.template getLastProp<VELOCITY_PREV>()[0] = 0.0;
		vectorDist.template getLastProp<VELOCITY_PREV>()[1] = 0.0;
		vectorDist.template getLastProp<VELOCITY_PREV>()[2] = 0.0;

		++obstacle_box;
	}

	vectorDist.map();

	// Now that we fill the vector with particles
	ModelCustom md;

	vectorDist.addComputationCosts(md);
	vectorDist.getDecomposition().decompose();
	vectorDist.map();

	///////////////////////////

	// Ok the initialization is done on CPU on GPU we are doing the main loop, so first we offload all properties on GPU

	vectorDist.hostToDevicePos();
	vectorDist.template hostToDeviceProp<TYPE,RHO,RHO_PREV,PRESSURE,VELOCITY,VELOCITY_PREV>();

	vectorDist.ghost_get<TYPE,RHO,PRESSURE,VELOCITY>(RUN_ON_DEVICE);

	auto cellList = vectorDist.getCellListGPU(2*H, CL_NON_SYMMETRIC, 2);

	timer tot_sim;
	tot_sim.start();

	size_t write = 0;
	size_t it = 0;
	size_t it_reb = 0;
	real_number t = 0.0;
	while (t <= t_end)
	{
		Vcluster<> & v_cl = create_vcluster();
		timer it_time;
		it_time.start();

		////// Do rebalancing every 200 timesteps
		it_reb++;
		if (it_reb == 300)
		{
			vectorDist.map(RUN_ON_DEVICE);

			// Rebalancer for now work on CPU , so move to CPU
			vectorDist.deviceToHostPos();
			vectorDist.template deviceToHostProp<TYPE>();

			it_reb = 0;
			ModelCustom md;
			vectorDist.addComputationCosts(md);
			vectorDist.getDecomposition().decompose();

			if (v_cl.getProcessUnitID() == 0)
			{std::cout << "REBALANCED " << it_reb << std::endl;}
		}

		vectorDist.map(RUN_ON_DEVICE);

		// Calculate pressure from the density
		EqState(vectorDist);

		real_number max_visc = 0.0;

		vectorDist.ghost_get<TYPE,RHO,PRESSURE,VELOCITY>(RUN_ON_DEVICE);

		// Calc forces
		calc_forces(vectorDist,cellList,max_visc,cnt);

		// Get the maximum viscosity term across processors
		v_cl.max(max_visc);
		v_cl.execute();

		// Calculate delta t integration
		real_number dt = calc_deltaT(vectorDist,max_visc);

		// VerletStep or euler step
		it++;
		if (it < 40)
			verlet_int(vectorDist,dt);
		else
		{
			euler_int(vectorDist,dt);
			it = 0;
		}

		t += dt;

		if (write < t*100)
		{
			// Sensor pressure require update ghost, so we ensure that particles are distributed correctly
			// and ghost are updated
			vectorDist.map(RUN_ON_DEVICE);
			vectorDist.ghost_get<TYPE,RHO,PRESSURE,VELOCITY>(RUN_ON_DEVICE);

			// calculate the pressure at the sensor points
			//sensor_pressure(vectorDist,cellList,press_t,probes);

			std::cout << "OUTPUT " << dt << std::endl;

			// When we write we have move all the particles information back to CPU

			vectorDist.deviceToHostPos();
			vectorDist.deviceToHostProp<TYPE,RHO,RHO_PREV,PRESSURE,DRHO,FORCE,VELOCITY,VELOCITY_PREV,RED,RED2>();

			vectorDist.write_frame("Geometry",write);
			write++;

			if (v_cl.getProcessUnitID() == 0)
			{std::cout << "TIME: " << t << "  write " << it_time.getwct() << "   " << it_reb << "   " << cnt << " Max visc: " << max_visc << "   " << vectorDist.size_local()  << std::endl;}
		}
		else
		{
			if (v_cl.getProcessUnitID() == 0)
			{std::cout << "TIME: " << t << "  " << it_time.getwct() << "   " << it_reb << "   " << cnt  << " Max visc: " << max_visc << "   " << vectorDist.size_local() << std::endl;}
		}
	}

	tot_sim.stop();
	std::cout << "Time to complete: " << tot_sim.getwct() << " seconds" << std::endl;

	openfpm_finalize();
}

#else

int main(int argc, char* argv[])
{
        return 0;
}

#endif
