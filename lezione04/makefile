LIBS2:= Vettore.o MolDyn_NVE.o random.o Experiment.o
LISTH:= MolDyn_NVE.h
CC = g++
CFLAGS = -Wall -O3 --std=c++11

%.o: %.cpp %.h ${LISTH}
	g++ -c $< ${CFLAGS}

main01.x: main01.C ${LIBS2}
	g++ $^ -o $@ ${CFLAGS}

main02.x: main02.C ${LIBS2}
	g++ $^ -o $@ ${CFLAGS}

main03.x: main03.C ${LIBS2}
	g++ $^ -o $@ ${CFLAGS}

clear:
	rm *.o *~ seed.out
