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
rnd->SetRandom("Primes","seed.in");	//parallel initialization
Population Pop(rnd,size,rank);

//simulation
Pop.Evolution();
cout<<endl;

MPI_Finalize();
return 0;
}
