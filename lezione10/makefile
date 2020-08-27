LIBS2:= random.o Vettore.o Experiment.o Metropolis.o TSPpath.o
CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

%.o: %.cpp %.h
	g++ -c $< ${CFLAGS}

GeneticTSP.o: GeneticTSP.cpp GeneticTSP.h
	${CC} -c $< ${CFLAGS}

main01.x: main01.C ${LIBS2} SimAnnTSP.o
	g++ $^ -o $@ ${CFLAGS}

main02.x: main02.C ${LIBS2} GeneticTSP.o
	${CC} $^ -o $@ ${CFLAGS}

clear:
	rm *.o *~
