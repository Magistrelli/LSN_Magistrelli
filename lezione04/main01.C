#include "MolDyn_NVE.h"

int main (int argc, char *argv[]){
if (argc!=3) {
    cerr << "Usage: " << argv[0] << " <Simulated state of matter> <restarting index (#restart-1)>" << endl;
    return -1;
}

//Variables
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");//initialization random_gen from files
MolDyn Test(1,rnd);		//simulation's initialization
const int nobs=4;		//#obs directed measured
string *state=new string(argv[1]);
string *iter= new string(argv[2]);
ofstream Results[nobs];

//compute
Results[0].open("results/output."+ *state +".epot."+ *iter);
Results[1].open("results/output."+ *state +".ekin."+ *iter);
Results[2].open("results/output."+ *state +".etot."+ *iter);
Results[3].open("results/output."+ *state +".temp."+ *iter);
Test.SetBurn(Test.GetL()*Test.GetBlk());
Test.InstantValues(Results,nobs,10,10);//run simulation, write results each 10 steps
Test.BoxScale();
Test.ConfFinalPlus();
for (int i=0; i<nobs; i++)	{Results[i].close();}

return 0;
}
