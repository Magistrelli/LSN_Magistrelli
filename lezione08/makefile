LIBS2:= random.o Vettore.o Experiment.o Metropolis.o VarMC.o
CC = g++
CFLAGS = -Wall -O3 --std=c++11

%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}

main01.x: main01.C ${LIBS2}
	$(CC) $^ -o $@ ${CFLAGS}

main02.x: main02.C ${LIBS2}
	$(CC) $^ -o $@ ${CFLAGS}


LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

qmc1d.o: qmc1d.cpp
	$(CC) -Wall -c $< ${INCS}

qmc1d.x: qmc1d.o
	$(CC) ${CFLAGS} -o $@ $^ ${LIBS}

qmc1d_const.o: qmc1d_const.cpp
	$(CC) -Wall -c $< ${INCS}

qmc1d_const.x: qmc1d_const.o
	$(CC) ${CFLAGS} -o $@ $^ ${LIBS}



clear:
	rm *.o *~
