#include "GeneticTSP.h"

int main (int argc, char *argv[]){

//parallelization stuff
int size,rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&size);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);

//Variables
Random* rnd = new Random();
rnd->SetParallel(size,rank);
rnd->SetRandom("Primes","seed.in");	//initialization parallel random_gen
Population Pop(rnd,size,rank);

//simulation
Pop.Evolution();			//comment if mutation check
cout<<endl;

MPI_Finalize();
return 0;
}
