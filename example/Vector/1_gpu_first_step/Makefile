include ../../example.mk
include ../../common.mk

OBJ = main.o

gpu_fstep:

%.o: %.cu
	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

gpu_fstep: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: gpu_fstep

run: gpu_fstep
	mpirun --oversubscribe -np 2 ./gpu_fstep

.PHONY: clean all run

clean:
	rm -f *.o *~ core gpu_fstep

