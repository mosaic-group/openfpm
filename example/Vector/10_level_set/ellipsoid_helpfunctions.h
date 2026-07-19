#include <math.h>
// This is a collection of helpfunctions that were used to run convergence tests and benchmarks for the ellipse/ellipsoid
// in the particle closest point draft. Computing the theoretical closest points and distances from a given query point
// is done using the first four functions which were adopted from David Eberly "Distance from a Point to an Ellipse, an 
// Ellipsoid, or a Hyperellipsoid", 2013.
//
// Created by lschulze
double GetRoot (double r0, double z0, double z1, double g)
    {
	const int maxIter = 100;
        double n0 = r0*z0;
        double s0 = z1 - 1;
        double s1 = ( g < 0 ? 0 : sqrt(n0*n0+z1*z1) - 1 ) ;
        double s = 0;
        for ( int i = 0; i < maxIter; ++i ){
	    if (i == (maxIter - 1)) std::cout<<"distance point ellipse algorithm did not converge."<<std::endl;
            s = ( s0 + s1 ) / 2 ;
            if ( s == s0 || s == s1 ) {break; }
            double ratio0 = n0 /( s + r0 );
            double ratio1 = z1 /( s + 1 );
            g = ratio0*ratio0 + ratio1*ratio1 - 1 ;
            if (g > 0) {s0 = s;} else if (g < 0) {s1 = s ;} else {break ;}
        }
        return s;
    }

double GetRoot(double r0, double r1, double z0, double z1, double z2, double g)
{	
	const int maxIter = 100;
	double n0 = r0*z0;
	double n1 = r1*z1;
	double s0 = z2 - 1;
	double s1 = (g < 0 ? 0 : sqrt(n0*n0 + n1*n1 + z2*z2) - 1) ;
	double s = s;
	for(int i = 0 ; i < maxIter ; ++i )
	{
		if (i == (maxIter - 1)) std::cout<<"distance point ellipse algorithm did not converge."<<std::endl;
		s = ( s0 + s1 ) / 2 ;
		if ( s == s0 || s == s1 ) {break; }
		double ratio0 = n0 / ( s + r0 ); 
		double ratio1 = n1 / ( s + r1 );
		double ratio2 = z2 / ( s + 1 );
		g = ratio0*ratio0 + ratio1*ratio1 +ratio2*ratio2 - 1;
		if ( g > 0 ) { s0 = s ;} 
		else if ( g < 0 ) { s1 = s ; }
		else {break;}
	}
	return (s);
}

double DistancePointEllipse(double e0, double e1, double y0, double y1, double& x0, double& x1)
    {
        double distance;
        if ( y1 > 0){
            if ( y0 > 0){
                double z0 = y0 / e0;
                double z1 = y1 / e1;
                double g = z0*z0+z1*z1 - 1;
                if ( g != 0){
                    double r0 = (e0/e1)*(e0/e1);
                    double sbar = GetRoot(r0 , z0 , z1 , g);
                    x0 = r0 * y0 /( sbar + r0 );
                    x1 = y1 /( sbar + 1 );
                    distance = sqrt( (x0-y0)*(x0-y0) + (x1-y1)*(x1-y1) );
                    }else{
                        x0 = y0;
                        x1 = y1;
                        distance = 0;
                    }
                }
                else // y0 == 0
                    {x0 = 0 ; x1 = e1 ; distance = abs( y1 - e1 );}
        }else{ // y1 == 0
            double numer0 = e0*y0 , denom0 = e0*e0 - e1*e1;
            if ( numer0 < denom0 ){
                    double xde0 = numer0/denom0;
                    x0 = e0*xde0 ; x1 = e1*sqrt(1 - xde0*xde0 );
                    distance = sqrt( (x0-y0)*(x0-y0) + x1*x1 );
                }else{
                    x0 = e0;
                    x1 = 0;
                    distance = abs( y0 - e0 );
            }
        }
        return distance;
    }

double DistancePointEllipsoid(double e0, double e1, double e2, double y0, double y1, double y2, double& x0, double& x1, double& x2)
{
	double distance;
	if( y2 > 0 )
	{
		if( y1 > 0 )
		{
			if( y0 > 0 )
			{
				double z0 = y0 / e0;
				double z1 = y1 / e1;
				double z2 = y2 / e2;
				double g = z0*z0 + z1*z1 + z2*z2 - 1 ;
				if( g != 0 )
				{
					double r0 = (e0/e2)*(e0/e2);
					double r1 = (e1/e2)*(e1/e2);
					double sbar = GetRoot ( r0 , r1 , z0 , z1 , z2 , g );
					x0 = r0 *y0 / ( sbar + r0 );
					x1 = r1 *y1 / ( sbar + r1 );
					x2 = y2 / ( sbar + 1 );
					distance = sqrt( (x0 - y0)*(x0 - y0) + (x1 - y1)*(x1 - y1) + (x2 - y2)*(x2 - y2));
				}
				else
				{
					x0 = y0;
					x1 = y1;
					x2 = y2;
					distance = 0;
				}
			}
			else // y0 == 0
			{
				x0 = 0;
				distance = DistancePointEllipse( e1 , e2 , y1 , y2, x1, x2);
			}
		}
		else // y1 == 0
		{
			if( y0 > 0 )
			{
				x1 = 0;
				distance = DistancePointEllipse( e0 , e2 , y0 , y2, x0, x2);
			}
			else // y0 == 0
			{
				x0 = 0;
				x1 = 0;
				x2 = e2;
				distance = abs(y2 - e2);
			}
		}
	}
	else // y2 == 0
	{
		double denom0 = e0*e0 - e2*e2;
		double denom1 = e1*e1 - e2*e2;
		double numer0 = e0*y0;
		double numer1 = e1*y1;
		bool computed = false;
		if((numer0 < denom0) && (numer1 < denom1))
		{
			double xde0 = numer0/denom0;
			double xde1 = numer1/denom1 ;
			double xde0sqr = xde0 *xde0;
			double xde1sqr = xde1 * xde1 ;
			double discr = 1 - xde0sqr - xde1sqr;
			if( discr > 0 )
			{
				x0 = e0*xde0;
				x1 = e1*xde1;
				x2 = e2*sqrt(discr);
				distance = sqrt((x0 - y0)*(x0 - y0) + (x1 - y1)*(x1 - y1) + x2*x2);
				computed = true;
			}
		}
		if( !computed )
		{
			x2 = 0;
			distance = DistancePointEllipse(e0 , e1 , y0 , y1, x0, x1);
		}
	}
	return distance;
}

constexpr unsigned int factorial(unsigned int x)
{
    unsigned int fact = 1;
    for(int i = 1; i < (x + 1); i++) fact = fact*i;
    return(fact);
}
constexpr unsigned int minter_lp_degree_one_num_coeffs(unsigned int dims, unsigned int poly_degree)
{
    return(factorial(dims + poly_degree)/(factorial(dims)*factorial(poly_degree)));
}

int return_sign(double phi)
{
	if (phi > 0) return 1;
	if (phi < 0) return -1;
	return 0;
}

double randMinusOneToOne()
{
    double temp = rand() / (RAND_MAX + 1.0);
    //std::cout<<(2.0*temp - 1.0)<<std::endl;
    return(2.0*temp - 1.0);
}

template <typename particles_type> inline void perturb_pos(particles_type & vd, double dp, double factor)
{
	const int dim = particles_type::dims;
	auto part = vd.getDomainIterator();

	while(part.isNext())
	{
		auto a = part.get();

		for (int k = 0; k<dim; k++)
		{
			vd.getPos(a)[k] += factor*randMinusOneToOne()*dp;
			//vd.getPos(a)[k] += factor*vd.template getProp<0>(a)*dp;
			//vd.getPos(a)[k] += factor*dp;
		}

		++part;
	}
}

template <typename particles_type, int sdf, int sdf_analytical, int cp_theoretical> inline void update_sdfs(particles_type & vd, double A, double B, double C)
{
	unsigned constexpr dim = particles_type::dims;
	auto part = vd.getDomainIterator();
	while (part.isNext())
	{
		auto a = part.get();
		Point<dim, double> xa = vd.getPos(a);
		if (dim == 2) 
		{
			vd.template getProp<sdf>(a) = - 1.0 + sqrt((xa[0]/A)*(xa[0]/A) + (xa[1]/B)*(xa[1]/B));
			//vd.template getProp<sdf>(a) = 10.0;
			double xcp_analytical;
			double ycp_analytical;
			vd.template getProp<sdf_analytical>(a) = return_sign(vd.template getProp<sdf>(a))*DistancePointEllipse(A, B, abs(xa[0]), abs(xa[1]), xcp_analytical, ycp_analytical);
			//debugging
			//vd.template getProp<sdf_analytical>(a) = 1.0 - sqrt((xa[0]/A)*(xa[0]/A) + (xa[1]/B)*(xa[1]/B));
			vd.template getProp<cp_theoretical>(a)[0] = return_sign(xa[0])*xcp_analytical;
			vd.template getProp<cp_theoretical>(a)[1] = return_sign(xa[1])*ycp_analytical;
		}
		else if (dim == 3) 
		{
			vd.template getProp<sdf>(a) = 1.0 - sqrt(pow(xa[0]/A, 2.0) + pow(xa[1]/B, 2.0) + pow(xa[2]/C, 2.0));
			double xcp_analytical;
			double ycp_analytical;
			double zcp_analytical;
			vd.template getProp<sdf_analytical>(a) = return_sign(vd.template getProp<sdf>(a))*DistancePointEllipsoid(A, B, C, abs(xa[0]), abs(xa[1]), abs(xa[2]), xcp_analytical, ycp_analytical, zcp_analytical);
			vd.template getProp<cp_theoretical>(a)[0] = return_sign(xa[0])*xcp_analytical;
			vd.template getProp<cp_theoretical>(a)[1] = return_sign(xa[1])*ycp_analytical;
			vd.template getProp<cp_theoretical>(a)[2] = return_sign(xa[2])*zcp_analytical;
		}
		// binary distribution
		//if ((- 1 + sqrt((x/A)*(x/A) + (y/B)*(y/B))) > 0) vd.getProp<sdf>(a) = dp;
		//else vd.getProp<sdf>(a) = -dp;

		// Saye's initial configuration
		//vd.getProp<sdf>(a) = (1 - exp(-(x-0.3)*(x-0.3)-(y-0.3)*(y-0.3)))*(sqrt(4*x*x+9*y*y)-1);
		++part;
	}
}

template <typename discretization_type, int sdf, int sdf_analytical> inline void get_interpol_error(discretization_type & vd, double narrow_band_half_width, double A, double B, double C)
{
	static constexpr int dim = discretization_type::dims;
	auto part = vd.getDomainIterator();
	double err = 0.0;
	double maxerr = 0.0;
	double cumerror = 0.0;
	int numparts = 0;
	while(part.isNext())
	{
		auto a = part.get();
		Point<dim, double> xa = vd.getPos(a);
		if (abs(vd.template getProp<sdf_analytical>(a)) > narrow_band_half_width)
		{
			++part;
			continue;
		}
		if (dim == 2) err = abs(vd.template getProp<sdf>(a) - (1.0 - sqrt((xa[0]/A)*(xa[0]/A) + (xa[1]/B)*(xa[1]/B))));
		else if (dim == 3) err = abs(vd.template getProp<sdf>(a) - (1.0 - sqrt((xa[0]/A)*(xa[0]/A) + (xa[1]/B)*(xa[1]/B) + (xa[2]/C)*(xa[2]/C))));
		else break;
		if (err > maxerr) maxerr = err;
		cumerror += err;
		numparts++;
		++part;
	}
	//std::cout<<"Average interpolation error is "<<cumerror/numparts<<"\nMax interpolation error is "<<maxerr<<std::endl;
	printf("%f, ",maxerr);

}

template <typename discretization_type, int sdf, int sdfgrad, int curvature, int cp, int sdf_analytical, int cp_theoretical> inline void get_max_error(discretization_type & vd, double narrow_band_half_width, double A, double B, double C)
{
	static constexpr int dim = discretization_type::dims;
	auto part = vd.getDomainIterator();
	double err;
	double maxerr = 0.0;
	double errkappa;
	double maxerrkappa = 0.0;
	double errnorm;
	double maxerrnorm = 0.0;
	double errcp;
	double maxerrcp = 0.0;
	double errone;
	double maxerrone = 0.0;
	double kappa = 0.0;
	double avgsdferr = 0.0;
	int county = 0;

	while(part.isNext())
	{
		auto a = part.get();
		// L infinity norm:
		//errnorm = abs(nxa[0]);
		//if (abs(nxa[1]) > errnorm) errnorm = abs(nxa[1]);
		//if (abs(nxa[2]) > errnorm) errnorm = abs(nxa[2]);
		if (abs(vd.template getProp<sdf_analytical>(a)) > narrow_band_half_width)
		{
			++part;
			continue;
		}
		Point<dim, double> norm_analytical;
		Point<dim, double> norm_computed;
		for (int k = 0; k < dim; k++) norm_computed[k] = vd.template getProp<sdfgrad>(a)[k];
		Point<dim, double> cp_computed;
		Point<dim, double> cp_theo;
		for (int k = 0; k < dim; k++)
		{
			cp_computed[k] = vd.template getProp<cp>(a)[k];
			cp_theo[k] = vd.template getProp<cp_theoretical>(a)[k];
		}
		err = abs(vd.template getProp<sdf>(a) - vd.template getProp<sdf_analytical>(a));
		if (vd.template getProp<4>(a))
		{
			avgsdferr += err;
			++county;
		}
		errcp = norm(cp_computed - cp_theo);
		
		if (dim == 2)
		{
			norm_analytical[0] = 2*cp_theo[0]/(A*A);
			norm_analytical[1] = 2*cp_theo[1]/(B*B);
			double norm_norm_analytical = norm(norm_analytical);
			for(int k = 0; k < dim; k++) norm_analytical[k] = return_sign(vd.template getProp<sdf>(a))*norm_analytical[k]/norm_norm_analytical;
			errnorm = norm(norm_analytical - norm_computed);

			kappa = A*B/std::pow(B/A*B/A*cp_theo[0]*cp_theo[0] + A/B*A/B*cp_theo[1]*cp_theo[1], 1.5);
		}
		else if (dim == 3)
		{
			norm_analytical[0] = -2*cp_theo[0]/(A*A);
			norm_analytical[1] = -2*cp_theo[1]/(B*B);
			norm_analytical[2] = -2*cp_theo[2]/(C*C);
			double norm_norm_analytical = norm(norm_analytical);
			for(int k = 0; k < dim; k++) norm_analytical[k] = return_sign(vd.template getProp<sdf>(a))*norm_analytical[k]/norm_norm_analytical;
			//std::cout<<norm(norm_analytical)<<std::endl;
			errnorm = norm(norm_analytical - norm_computed);

			kappa = -(std::abs(cp_theo[0]*cp_theo[0] + cp_theo[1]*cp_theo[1] + cp_theo[2]*cp_theo[2] - A*A - B*B - C*C))/(2*A*A*B*B*C*C*std::pow(cp_theo[0]*cp_theo[0]/std::pow(A, 4) + cp_theo[1]*cp_theo[1]/std::pow(B, 4) + cp_theo[2]*cp_theo[2]/std::pow(C, 4), 1.5));
		}
		errkappa = abs(vd.template getProp<curvature>(a) - kappa);
		if (err > maxerr) maxerr = err;
		if (errkappa > maxerrkappa) maxerrkappa = errkappa;
		if (errnorm > maxerrnorm) maxerrnorm = errnorm;
		if (errcp > maxerrcp) maxerrcp = errcp;
		++part;
	}
	std::cout<<"Average error for sdf in narrow band is: "<<avgsdferr/county<<std::endl;
	std::cout<<"Maximum error for sdf on processor is: "<<maxerr<<std::endl;
	std::cout<<"Maximum error for surface normal on processor is: "<<maxerrnorm<<std::endl;
	std::cout<<"Maximum error for curvature on processor is: "<<maxerrkappa<<std::endl;
	std::cout<<"Maximum error for cp(xa) is: "<<maxerrcp<<std::endl;
	//printf("%.20f,\n ", maxerr);
}
