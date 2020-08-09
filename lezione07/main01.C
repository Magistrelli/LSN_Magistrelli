#include "Monte_Carlo_NVT.h"

int main (int argc, char *argv[]){
if (argc!=2) {
    cerr << "Usage: " << argv[0] << " <Simulated state of matter>" << endl;
    return -1;
}

//Variables
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");//initialization random_gen from files
CanonicEns Test(rnd);
const int nobs=2;		//obs:V,P (for g not blocking graph)
string *nameFile,*state=new string(argv[1]);
ofstream OutBurn[nobs],OutInst[nobs],OutFinal;

//equilibration
OutBurn[0].open("results/burnin."+ *state +".epot.ist");
OutBurn[1].open("results/burnin."+ *state +".pres.ist");
Test.SetBurn(int(1e4));
Test.InstantValues(OutBurn,nobs,1,100);//print also XYZ config every 100 steps
for(int i=0; i<nobs; i++) {OutBurn[i].close();}

//simulation
OutInst[0].open("results/output."+ *state +".epot.ist");
OutInst[1].open("results/output."+ *state +".pres.ist");
Test.SetBurn(int(5e5));
Test.InstantValues(OutInst,nobs,1,0);//obs value every step, not XYZ config
for(int i=0; i<nobs; i++) {OutInst[i].close();}

nameFile=new string("finalConf/config."+*state+".final");
cout << "Print final configuration to file " << *nameFile << endl << endl;
OutFinal.open(*nameFile);
Test.BoxScale();		//position conversion in Box lenght units
Test.ConfFinal(OutFinal);
OutFinal.close();

return 0;
}
