include ../../example.mk
include ../../common.mk

OBJ_DORD = main.o

all: complex_use

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

complex_use: $(OBJ_DORD)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all
	mpirun --oversubscribe -np 3 ./complex_use

.PHONY: clean all run all_test on_test

clean:
	rm -f *.o *~ core complex_use

