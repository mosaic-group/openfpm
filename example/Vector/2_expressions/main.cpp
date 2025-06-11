#include "Vector/vector_dist.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"

/*!
 * \page Vector_2_expression Vector 2 expression
 *
 * [TOC]
 *
 * # Vector expressions # {#e2_ve}
 * 
 * This example show how to use expressions to do complex calculation on particles
 * 
 * ## Initialization ## {#e2_ve_init}
 *
 * He define some useful constants
 *
 * \snippet Vector/2_expressions/main.cpp Initialization and parameters
 *
 * 
 */

//! \cond [Initialization and parameters] \endcond

constexpr int A = 0;
constexpr int B = 1;
constexpr int C = 2;
constexpr int D = 3;

//! \cond [Initialization and parameters] \endcond


/*!
 * \page Vector_2_expression Vector 2 expressions
 *
 * ## Kernel function ## {#e2_ve_ker}
 *
 * We define a kernel function. In general a kernel is a function like \f$ K(x_p,x_q,A_p,A_q) \f$.
 * Where \f$ x_p \f$ and \f$ x_q \f$ are positions of the 2 particles p and q. \f$ A_p \f$ and \f$ A_q \f$ are
 * the properties on p and q.
 *
 * A second version more General is
 *
 * \f$ K(p,q,A_p,A_q,particles) \f$
 *
 * Where \f$ p \f$ in the index of the particle and \f$ q \f$ are the index of the particles produced by the cell-List.
 * and particles is the vector containing the particle p and q (In multi-phase SPH this can also be not true)
 *
 * \note On multi-phase this is not true
 *
 * In this case we are defining an exponential-like function \f$ A_q e^{\frac{|x_p-x_q|^2}{\sigma}} \f$
 *
 * \snippet Vector/2_expressions/main.cpp exp ker
 *
 *
 */

//! \cond [exp ker] \endcond

// Kernel are structures
struct exp_kernel
{
	// variance of the exponential
	float sigma;

	// Constructor require to define the variance
	exp_kernel(float sigma)
	:sigma(sigma)
	{}

	// The kernel function itself K(x_p,x_q,A_p,A_q). The name MUST be value and require 4 arguments.
	// p,q,A_p,A_q Position of p, position of q, value of the property A on p, value of the property A on q
	// and it return the same type of the property A
	inline Point<3,double> value(const Point<3,float> & p, const Point<3,float> & q,Point<3,double> & Pp,Point<3,double> & Pq)
	{
		// calculate the distance between p and q
		float dist = norm(p-q);

		// Calculate the exponential
		return Pq * exp(- dist * dist / sigma);
	}

	// The kernel function itself K(p,q,A_p,A_q,particles). The name MUST be value and require 5 arguments.
	// p,q,A_p,A_q,particles size_t index of p, size_t index of q, value of the property A on p, value of the property A on q,
	// original vector of the particles p, and it return the same type of the property A
	inline Point<3,double> value(size_t p, size_t q, Point<3,double> & Pp, Point<3,double> & Pq, const vector_dist<3,double, aggregate<double,double,Point<3,double>,Point<3,double>> > & v)
	{
		// calculate the distance between p and q
		float dist = norm(p-q);

		// Calculate the exponential
		return Pq * exp(- dist * dist / sigma);
	}
};

//! \cond [exp ker] \endcond

int main(int argc, char* argv[])
{
	/*!
	 * \page Vector_2_expression Vector 2 expressions
	 *
	 * ## Initialization ## {#e2_ex_init}
	 *
	 * Here we Initialize the library, we create a Box that define our domain, boundary conditions, and ghost
	 *
	 * \see \ref e0_s_init
	 *
	 * \snippet Vector/2_expressions/main.cpp init
	 *
	 */
	
	//! \cond [init] \endcond

	// initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with intervals from 0 to 1.0
	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Here we define the boundary conditions of our problem
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<3,double> g(0.01);

	// Delta t
	double dt = 0.01;

	//! \cond [init] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expressions
	 *
	 * ## %Vector create ##
	 *
	 * We create 2 vectors one with 4 properties 2 scalar and 2 vectors. and one with one scalar
	 *
	 * \see \ref e0_s_vector_inst
	 *
	 * \snippet Vector/2_expressions/main.cpp inst
	 *
	 */

	//! \cond [inst] \endcond

	// distributed vectors
	vector_dist<3,double, aggregate<double,double,Point<3,double>,Point<3,double>> > vd(4096,domain,bc,g);
	vector_dist<3,double, aggregate<double> > vd2(4096,domain,bc,g);

	//! \cond [inst] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expressions
	 *
	 * ## Loop over particles ##
	 *
	 * We assign random position to particles for both vectors
	 *
	 * \see \ref e0_s_assign_pos
	 *
	 * \snippet Vector/2_expressions/main.cpp assign
	 *
	 */

	//! \cond [assign] \endcond

	// Assign random position to the vector dist
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.getPos(p)[0]  = (double)rand() / RAND_MAX;
		vd.getPos(p)[1]  = (double)rand() / RAND_MAX;
		vd.getPos(p)[2]  = (double)rand() / RAND_MAX;

		vd2.getPos(p)[0] = (double)rand() / RAND_MAX;
		vd2.getPos(p)[1]  = (double)rand() / RAND_MAX;
		vd2.getPos(p)[2]  = (double)rand() / RAND_MAX;

		++it;
	}

	//! \cond [assign] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expression
	 *
	 * ## Mapping particles ## {#e2_ex_map}
	 *
	 * Redistribute the particles according to the decomposition
	 *
	 * \see \ref e0_s_map
	 *
	 * \snippet Vector/2_expressions/main.cpp map
	 *
	 */

	//! \cond [map] \endcond

	vd.map();
	vd.map();

	//! \cond [map] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expression
	 *
	 * ## Getting expression properties ## {#e2_ex_epp}
	 *
	 * Before using expression we have to get alias for each properties we want to use for the expression
	 *
	 *
	 * \snippet Vector/2_expressions/main.cpp expr init
	 *
	 */

	//! \cond [expr init] \endcond

	// vA is an alias for the property A of the vector vd
	// vB is an alias for the property V of the vector vd
	// vC is ...
	auto vA = getV<A>(vd);
	auto vB = getV<B>(vd);
	auto vC = getV<C>(vd);
	auto vD = getV<D>(vd);

	// This indicate the particle position
	auto vPOS = getV<POS_PROP>(vd);

	// same concept for the vector v2
	auto v2A = getV<0>(vd);
	auto v2POS = getV<POS_PROP>(vd);

	//! \cond [expr init] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expression
	 *
	 * ## Expressions ## {#e2_ex_use}
	 *
	 * Once we have the expression properties we can use them to compose expressions
	 *
	 *
	 * \snippet Vector/2_expressions/main.cpp example expr
	 *
	 */

	//! \cond [example expr] \endcond

	// Assign 1 to the property A
	vA = 1;

	// Equal the property A to B
	vB = vA;

	// Assign to the property C and for each component 4.0
	vC = 4;

	// Assign the position of the particle to the property D
	vD = vPOS;

	//! \cond [example expr] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expression
	 *
	 *
	 * All these expression are applied point wise for each particles
	 *
	 *
	 *
	 * \note the * is a scalar product in case of vectors
	 * 		 simple multiplication in case of scalar
	 *
	 *
	 * \snippet Vector/2_expressions/main.cpp pw expr
	 *
	 */

	//! \cond [pw expr] \endcond

	// In vA we set the distance from the origin for each particle
	vA = norm(vPOS);

	//
	// For each particle p calculate the expression under
	//
	// NOTE sin(2.0 * vD)  and exp(5.0 * vD) are both vector
	//
	//              so this is a scalar product
	//                      |
	//                      V
	vA = vA + sin(2.0 * vD) * exp(1.2 * vD);
    //       |_____________________________|
    //                      |
    //                      V
    //              and this is a number
    //
    //


	// pmul indicate component-wise multiplication the same as .* in Matlab
	vC = pmul(vD,vD) * dt;

	// Normalization of the vector vD
	vD = vD / sqrt( vD * vD );

	//! \cond [pw expr] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expression
	 *
	 * ## Visualization, write VTK files ##
	 *
	 * With this function we output the particles in VTK format.
	 *
	 * \see \ref e0_s_vis_vtk
	 *
	 * \snippet Vector/2_expressions/main.cpp vtk
	 *
	 */

	//! \cond [vtk] \endcond

	// we write the result on file
	vd.write("output");

	//! \cond [vtk] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expression
	 *
	 * ## Lazy expression and convert object into expressions ##
	 *
	 * If we want to use an object like Point into an expression we have to first convert into an expression
	 * with the function getVExpr. Every expression does not produce any code execution until we do not assin
	 * this expression to some property
	 *
	 * \snippet Vector/2_expressions/main.cpp lazy expr
	 *
	 */

	//! \cond [lazy expr] \endcond

	// Constant point
	Point<3,double> p0({0.5,0.5,0.5});

	// we cannot use p0 directly we have to create an expression from it
	auto p0_e = getVExpr(p0);

	// here we are creating an expression, expr1 collect the full expression but does not produce any
	// code execution
	auto expr1 = 1.0/2.0/sqrt(M_PI)*exp(-(v2POS-p0_e)*(v2POS-p0_e)/2.0);

	// here the expression is executed on all particles of v2 and the result assigned to the property A of v2
	v2A = expr1;

	//! \cond [lazy expr] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expression
	 *
	 *
	 * Each operator= produce an iteration across particles. It is possible to use the function assign
	 * to execute multiple assignment in one cycle. Assign support a variable number of expression
	 *
	 * \snippet Vector/2_expressions/main.cpp assign expr
	 *
	 */

	//! \cond [assign expr] \endcond

	assign(vA,1.0,               // vA = 1
		   vB,expr1,          // vB = 1.0/2.0/sqrt(M_PI)*exp(-(v2POS-p0_e)*(v2POS-p0_e)/2.0)
		   vC,vD/norm(vD),    // vC = vD/norm(vD)
		   vD,2.0*vC);        // vD = 2.0*vC // here vC = vD/norm(vD)

	//! \cond [assign expr] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expression
	 *
	 * ## Apply kernel ##
	 *
	 *
	 * we shown the simple usage of simple point-wise function like sin
	 * and exp. Here we show the function applyKernel that encapsulate
	 * the operation
	 *
	 * \f$\sum_{q = Neighborhood(p)} Kernel(x_p,x_q,A_p,A_q) \f$
	 *
	 * where \f$ x_p \f$ is the position of the particle p,
	 *       \f$ x_q \f$ is the position of the particle q,
	 *       \f$ A_p \f$ is the value of the property carried by the particle p,
	 *       \f$ A_q \f$ is the value of the property carried by the particle q
	 *
	 * \note On particle method we saw the following formulas
	 *
	 * \note \f$  \sum_{q = Neighborhood(p)} A_q [D^{\beta}ker](x_p,x_q) V_q \f$ in SPH
	 *
	 * \note \f$  \sum_{q = Neighborhood(p)} (A_p-A_q)ker(x_p,x_q) V_q \f$ in PSE
	 *
	 * \note These 2 formulas and others can be embedded into the previous one where
	 *
	 * \note \f$ Kernel(x_p,x_q,A_p,A_q) = A_q [D^{\beta}ker](x_p,x_q) V_q \f$
	 *
	 * \note \f$ Kernel(x_p,x_q,A_p,A_q) = (A_p-A_q)*ker(x_p,x_q) V_q \f$
	 *
	 * We create some useful object to use apply kernel plus we synchronize the ghost for the
	 * property 3 because apply kernel require to read the ghost part
	 *
	 * \snippet Vector/2_expressions/main.cpp apply ker expr init
	 *
	 */

	//! \cond [apply ker expr init] \endcond

	// Create a kernel with sigma 0.5
	exp_kernel ker(0.5);

	// Get a Cell list with 0.1 of cutoff-radius
	auto cl = vd.getCellList(0.1);

	// We are going to do some calculation that require ghost, so sync it 3 == vD
    vd.ghost_get<3>();

    //! \cond [apply ker expr init] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expression
     *
     * Than we construct an expression with applyKernel
     *
     * \snippet Vector/2_expressions/main.cpp apply ker expr
     *
	 * A renormalizaion of the vector vD is calculated on fly for each particle q neighborhood of p. The
	 * the following expression is calculated for each particle p
	 *
	 * \f$ (\sum_{q = Neighborhood p} Ker(x_p,x_q,(\frac{vD}{norm(vD)})_p,(\frac{vD}{norm(vD)})_q)) + vD \f$
	 *
	 *
	 * \note Cerefull here ApplyKernel is a special function if it act on property vD we cannot write on property vD
	 *  ( vD = ...  is an error and produce unpredictable result )
	 *
	 * Exist also a second version suitable for more General cases
	 *
	 * \snippet Vector/2_expressions/main.cpp apply ker expr gen
	 *
	 * In the case of General-Finite-Differences, DC-PSE and Multi-Phase SPH. \f$ K(x_p,x_q,A_p,A_q)\f$ is limitative
	 * and we have to use the following \f$ K(p,q,A_p,A_q,particles) \f$.
	 * In this case the position \f$ x_p \f$ and \f$ x_q \f$ are substituted by the id of p and id of q comming from
	 * the cell list + the vector on which we are operating.
	 *
	 * \note in case of Multi-Phase multiple vectors are needed. But such vector can passed in the Kernel structure
	 *
	 * \snippet Vector/2_expressions/main.cpp apply ker expr gen
	 *
	 */

    //! \cond [apply ker expr] \endcond

	vC = applyKernel_in( vD / norm(vD) ,vd,cl,ker) + vD;

	//! \cond [apply ker expr] \endcond

	//! \cond [apply ker expr gen] \endcond

	// Second version in this case the kernel is called with the particles id for p
	// and
	vC = applyKernel_in_gen( vD / norm(vD) ,vd,cl,ker) + vD;

	//! \cond [apply ker expr gen] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expression
     *
     * If the expression is point-wise the rule the propety on the right cannot appear on the left does not
     * apply
     *
     */

	// ok, because does not use applyKernel
	vD = vD / norm(vD);

	//! \cond [apply ker expr] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expression
	 *
	 * ## VTK and visualization ##
	 *
	 * Write a VTK file and visualize. Before write the particles also delete Ghost
	 * (VTK writer in general output also the ghost particles)
	 *
	 * \see \ref e0_s_vis_vtk
	 *
	 *
	 * \snippet Vector/1_celllist/main.cpp finalize
	 *
	 */

	//! \cond [vtk] \endcond

	vd.deleteGhost();
	vd.write("output2");

	//! \cond [vtk] \endcond

	/*!
	 * \page Vector_2_expression Vector 2 expression
	 *
	 * ## Finalize ## {#finalize}
	 *
	 *  At the very end of the program we have always de-initialize the library
	 *
	 * \snippet Vector/0_simple/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond
}
