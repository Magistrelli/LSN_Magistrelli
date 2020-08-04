#include "MolDyn_NVE.h"

int main (int argc, char *argv[]){
if (argc!=2) {
    cerr << "Usage: " << argv[0] << " <Simulated state of matter>" << endl;
    return -1;
}

//Variables
int nblk=int(1e2),nobs=5;	//number of blocks and of obs
int nRestart=8;			//number of restarts
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");//initialization random_gen from files
MolDyn Test(nblk,rnd);		//simulation's initialization
string *state=new string(argv[1]);
ofstream Results[nobs],Null;

//compute
Test.SetBurn(int(1e3));
for(int i=0; i<nRestart; i++){	//warm-up, convergence to equilibrium
    cout << "Ignition number: " << i+1 << endl;
    Test.BurnIn(Null,0);	//no output
    Test.Restart();
}
cout << endl;

Results[0].open("results/Ar.ave."+ *state +".epot.out");
Results[1].open("results/Ar.ave."+ *state +".ekin.out");
Results[2].open("results/Ar.ave."+ *state +".etot.out");
Results[3].open("results/Ar.ave."+ *state +".temp.out");
Results[4].open("results/Ar.ave."+ *state +".gofr.out");
Test.DoSampling(Results,Null,0);//run simulation
for (int k=0; k<nobs; k++)	{Results[k].close();}

return 0;
}
