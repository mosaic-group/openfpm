include ../../example.mk
include ../../common.mk

OBJ_GPU = main_gpu.o
OBJ_CPU = main_cpu.o
OBJ_CPU_BEST = main_cpu_best.o

all: md_dyn_gpu md_dyn_cpu md_dyn_cpu_best

md_dyn_gpu_test: OPT += -DTEST_RUN
md_dyn_gpu_test: all

%.o: %.cu
	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

md_dyn_gpu: $(OBJ_GPU)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

md_dyn_cpu: $(OBJ_CPU)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

md_dyn_cpu_best: $(OBJ_CPU_BEST)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: md_dyn_gpu_test
	mpirun --oversubscribe -np 3 ./md_dyn_gpu && mpirun --oversubscribe -np 3 ./md_dyn_cpu && mpirun --oversubscribe -np 3 ./md_dyn_cpu_best;

.PHONY: clean all run

clean:
	rm -f *.o *~ core md_dyn_gpu md_dyn_cpu md_dyn_cpu_best

