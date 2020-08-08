LIBS2:= random.o Vettore.o Experiment.o
LIBSM:= Metropolis.o Monte_Carlo_NVT.o
CC = g++
CFLAGS = -Wall -O3 --std=c++11

%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}

main01.x: main01.C ${LIBS2} ${LIBSM}
	$(CC) $^ -o $@ ${CFLAGS}

main02.x: main02.C ${LIBS2} ${LIBSM}
	$(CC) $^ -o $@ ${CFLAGS}

main03.x: main03.C ${LIBS2} MolDyn_NVE.o
	$(CC) $^ -o $@ ${CFLAGS}

clear:
	rm *.o *~
