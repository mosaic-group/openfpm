include ../../example.mk
include ../../common.mk

OBJ_VIC_PETSC = main_vic_petsc.o
OBJ_VIC_PETSC_OPT = main_vic_petsc_opt.o

all: vic_petsc vic_petsc_opt

vic_petsc_test: OPT += -DTEST_RUN
vic_petsc_test: vic_petsc vic_petsc_opt

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

vic_petsc_opt: $(OBJ_VIC_PETSC_OPT)
	$(CC) -o $@ $^ $(LIBS_PATH) $(LIBS)

vic_petsc: $(OBJ_VIC_PETSC)
	$(CC) -o $@ $^ $(LIBS_PATH) $(LIBS)
	
run: vic_petsc_test
	mpirun --oversubscribe -np 4 ./vic_petsc && mpirun --oversubscribe -np 4 ./vic_petsc_opt
	
.PHONY: clean all

clean:
	rm -f *.o *~ core vic

