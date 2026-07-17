
#ifdef __NVCC__

#define PRINT_STACKTRACE
//#define STOP_ON_ERROR
#define OPENMPI
//#define SE_CLASS1

//#define USE_LOW_REGISTER_ITERATOR

#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include "util/stat/common_statistics.hpp"


__global__ void test1_syncthreads()
{
    __syncthreads();
    __syncthreads();
    __syncthreads();
    __syncthreads();
    __syncthreads();
    __syncthreads();

    __syncthreads();
    __syncthreads();
    __syncthreads();
    __syncthreads();
    __syncthreads();
    __syncthreads();

    __syncthreads();
    __syncthreads();
    __syncthreads();
    __syncthreads();
    __syncthreads();
    __syncthreads();

    __syncthreads();
    __syncthreads();
    __syncthreads();
    __syncthreads();
    __syncthreads();
    __syncthreads();
}


struct ite_g
{
    dim3 wthr;
    dim3 thr;

    size_t nblocks() const
	{
		return wthr.x * wthr.y * wthr.z;
	}

	size_t nthrs() const
	{
		return thr.x * thr.y * thr.z;
	}
};

int main(int argc, char* argv[])
{

    // initialize the library
	openfpm_init(&argc,&argv);

	openfpm::vector<double> tele_ker;

	ite_g g;

    	g.wthr = dim3(512*512,1,1);
    	g.thr = dim3(8,1,1);

	for (int i = 0; i < 10; i++)
	{
            timer t_ker;
            t_ker.start();

            CUDA_LAUNCH(test1_syncthreads,g);

            t_ker.stop();

	    std::cout << "TKERNEL: " << t_ker.getwct() << std::endl;
	    



	//////////////////////////////////////////////////////////////////////////////////////////////////

	    tele_ker.add(t_ker.getwct());


      ///////////////////////////////////////////////////////////////////////////////////////////////////////

	}

	double tele_ker_mean;
        double tele_ker_dev;
	standard_deviation(tele_ker,tele_ker_mean,tele_ker_dev);

	std::cout << g.wthr.x*g.wthr.y*g.wthr.z << "  " << g.thr.x << std::endl;
	std::cout << "SYNCTHREAD LATENCY: " << tele_ker_mean / (g.wthr.x*g.wthr.y*g.wthr.z*24*g.thr.x) * 1e9 << " ns " << " error: " << tele_ker_dev << std::endl;

	openfpm_finalize();
}
 
#else

int main(int argc, char* argv[])
{
        return 0;
}

#endif
