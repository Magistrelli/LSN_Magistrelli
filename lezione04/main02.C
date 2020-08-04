#include "MolDyn_NVE.h"

int main (int argc, char *argv[]){
if (argc!=2) {
    cerr << "Usage: " << argv[0] << " <Simulated state of matter>" << endl;
    return -1;
}

//Variables
int nblk=int(1e2);		//number of blocks
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");//initialization random_gen from files
MolDyn Test(nblk,rnd);		//simulation's initialization
string *state=new string(argv[1]);
const int nobs=5;
ofstream Results[nobs],Null;

//compute
Results[0].open("results/ave."+ *state +".epot.out");
Results[1].open("results/ave."+ *state +".ekin.out");
Results[2].open("results/ave."+ *state +".etot.out");
Results[3].open("results/ave."+ *state +".temp.out");
Results[4].open("results/ave."+ *state +".gofr.out");
Test.DoSampling(Results,Null,0);		//run simulation
for (int i=0; i<nobs; ++i)	{Results[i].close();}

return 0;
}
