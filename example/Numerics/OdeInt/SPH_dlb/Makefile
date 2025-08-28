include ../../../example.mk
include ../../../common.mk

OBJ = main.o
OBJ2 = main2.o

sph_dlb_test: OPT += -DTEST_RUN
sph_dlb_test: sph_dlb sph_dlb_odeint

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

sph_dlb: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

sph_dlb_odeint: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: sph_dlb sph_dlb_odeint

run: sph_dlb_test
	mpirun --oversubscribe -np 2 ./sph_dlb && mpirun --oversubscribe -np 2 ./sph_dlb_odeint

.PHONY: clean all run

clean:
	rm -f *.o *~ core sph_dlb sph_dlb_odeint 

