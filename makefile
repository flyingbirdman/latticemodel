#-----------------------------------
#configure
#-----------------------------------

OFLAGS = -O3 
LIBFLAGS = -llapacke -llapack -lrefblas
MC = g++

#-----------------------------------
#List objects
#-----------------------------------

OBJS = nb.o

#-----------------------------------
#Object compilation rules
#-----------------------------------
nb.o : nb.cpp nb.h nrlib.h vector_operation.h matrix_operation.h
	$(MC) $(OFLAGS) -c nb.cpp

2dlattice.o : 2dlattice.cpp 2dlattice.h nrlib.h nb.o
	$(MC) $(OFLAGS) -c 2dlattice.cpp

1dlattice.o : 1dlattice.cpp 1dlattice.h nrlib.h nb.o
	$(MC) $(OFLAGS) -c 1dlattice.cpp

test_2dlattice.o : test_2dlattice.cpp 2dlattice.o
	$(MC) $(OFLAGS) -c test_2dlattice.cpp

test_1dlattice.o : test_1dlattice.cpp 1dlattice.o
	$(MC) $(OFLAGS) -c test_1dlattice.cpp

test_d2_lattice.o : test_d2_lattice.cpp nb.o
	$(MC) $(OFLAGS) -c test_d2_lattice.cpp

test_d1_lattice.o : test_d1_lattice.cpp nb.o
	$(MC) $(OFLAGS) -c test_d1_lattice.cpp

#-----------------------------------
#Primary rules
#-----------------------------------

1dtest : test_d1_lattice.o $(OBJS)
	$(MC) $(OFLAGS) -o execute_test_d1 test_d1_lattice.o $(LIBFLAGS)

2dtest : test_d2_lattice.o $(OBJS)
	$(MC) $(OFLAGS) -o execute_test_d2 test_d2_lattice.o $(LIBFLAGS)

2dlattice : test_2dlattice.o $(OBJS)
	$(MC) $(OFLAGS) -o execute_test_2dlattice test_2dlattice.o $(LIBFLAGS)

1dlattice : test_1dlattice.o $(OBJS)
	$(MC) $(OFLAGS) -o execute_test_1dlattice test_1dlattice.o $(LIBFLAGS)


clean: #removes all old *.o and *.mod files
	rm *.o
	rm execute*
