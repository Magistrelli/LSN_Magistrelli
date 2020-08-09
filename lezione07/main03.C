#include "MolDyn_NVE.h"

int main (int argc, char *argv[]){
if (argc!=3) {
    cerr << "Usage: " << argv[0] << " <Simulated state of matter> <number of wanted restarts>" << endl;
    return -1;
}

//Variables
int nblk=20;			//#blk (as in input of NVT, 1000 stp/blk, #tot steps as NVT)
int nRestart=atoi(argv[2]);	//#restarts
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");//initialization random_gen from files
MolDyn Test(nblk,rnd);		//simulation's initialization
const int nbins=Test.GetNBins(), igr=Test.Get_iGr();
double* resGr = new double[2];		//couple of double (res,err)
string *state=new string(argv[1]),*nres=new string(argv[2]);
ofstream Results[igr+1],OutHisto,Null;
const int wd=20;

//compute
Test.SetBurn(int(1e3));
for(int i=0; i<nRestart; i++){	//warm-up, convergence to equilibrium
    cout << "Ignition number: " << i+1 << endl;
    Test.BurnIn(Null,0);	//no output
    Test.Restart();
}
cout << endl;

Results[igr].open("results/output."+*state+".gofr."+*nres+".NVE");
Test.DoSampling(Results,Null,0);//run simulation, write only g(r)
Results[igr].close();

OutHisto.open("results/output."+ *state +".gave."+*nres+".NVE");
OutHisto << Test.GetL() << setw(wd) << Test.GetBinWidth() << endl << endl;

for(int ibin=0; ibin<nbins; ibin++){
    resGr=Test.GetResErr(igr+ibin);
    OutHisto<<setw(wd)<< ibin <<setw(wd)<< resGr[0] <<setw(wd)<< resGr[1]<<endl;
}
OutHisto.close();

return 0;
}
