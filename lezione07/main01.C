#include "Monte_Carlo_NVT.h"

int main (int argc, char *argv[]){
if (argc!=2) {
    cerr << "Usage: " << argv[0] << " <Simulated state of matter>" << endl;
    return -1;
}

//Variables
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");	//initialization random_gen from files
CanonicEns Test(rnd);
const int nobs=2;			//obs:V,P (for g not blocking graph)
string *nameFile,*state=new string(argv[1]);
ofstream OutBurn[nobs],OutFinal;

//simulation
OutBurn[0].open("results/output."+ *state +".epot.ist");
OutBurn[1].open("results/output."+ *state +".pres.ist");
Test.SetBurn(int(5e5));
Test.InstantValues(OutBurn,nobs,1,0);
for(int i=0; i<nobs; i++) {OutBurn[i].close();}

nameFile=new string("finalConf/config."+*state+".final");
cout << "Print final configuration to file " << *nameFile << endl << endl;
OutFinal.open(*nameFile);
Test.BoxScale();
Test.ConfFinal(OutFinal);
OutFinal.close();

return 0;
}
