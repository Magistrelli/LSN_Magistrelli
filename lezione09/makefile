LIBS2:= random.o Vettore.o TSPpath.o GeneticTSP.o
CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

%.o: %.cpp %.h
	g++ -c $< ${CFLAGS}

GeneticTSP.o: GeneticTSP.cpp GeneticTSP.h
	${CC} -c $< ${CFLAGS}

main01.x: main01.C ${LIBS2}
	${CC} $^ -o $@ ${CFLAGS}

clear:
	rm *.o *~
