include ../../example.mk
include ../../common.mk

OBJ = main.o

all: hdf5


%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

hdf5: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)


run: hdf5
	mpirun --oversubscribe -np 2 ./hdf5

.PHONY: clean all run

clean:
	rm -f *.o *~ core hdf5

