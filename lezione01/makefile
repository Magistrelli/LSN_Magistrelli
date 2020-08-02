LIBS2:= random.o Vettore.o
CC = g++
CFLAGS = -Wall -O3 --std=c++11

%.o: %.cpp %.h
	g++ -c $< ${CFLAGS}

main01.x: main01.C ${LIBS2}
	g++ $^ -o $@ ${CFLAGS}

main02.x: main02.C ${LIBS2}
	g++ $^ -o $@ ${CFLAGS}

main03.x: main03.C ${LIBS2}
	g++ $^ -o $@ ${CFLAGS}

clear:
	rm *.o *~
