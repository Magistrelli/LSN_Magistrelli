#include "Monte_Carlo_NVT.h"

int main (int argc, char *argv[]){
if (argc!=3) {
    cerr << "Usage: " << argv[0] << " <Simulated state of matter> <restarting index (#restart-1)>" << endl;
    return -1;
}

//Variables
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");	//initialization random_gen from files
CanonicEns Test(rnd);
const int nobs=2,nbins=Test.GetNBins();	//obs: V, P, g(r)
double* resGr = new double[2];		//couple of double (res,err)
ofstream Null,OutRes[nobs+1],OutFinal,OutHisto;
string *state=new string(argv[1]),*nameFile;
string *iter= new string(argv[2]);
const int wd=20;

//simulation
Test.SetBurn(int(3e3));
Test.BurnIn(Null,0);

OutRes[0].open("results/output."+ *state +".epot."+ *iter);
OutRes[1].open("results/output."+ *state +".pres."+ *iter);
OutRes[2].open("results/output."+ *state +".gofr."+ *iter);
Test.DoSampling(OutRes,Null,0);
for(int i=0; i<nobs+1; i++) {OutRes[i].close();}

nameFile=new string("finalConf/config."+*state+".final."+*iter);
cout << "Print final configuration to file "+*nameFile << endl << endl;
OutFinal.open(*nameFile);
Test.BoxScale();			//pos conversion in Box lenght units
Test.ConfFinal(OutFinal);
OutFinal.close();

OutHisto.open("results/output."+ *state +".gave."+ *iter);
OutHisto << Test.GetL() << setw(wd) << Test.GetBinWidth() << endl << endl;
for(int ibin=0; ibin<nbins; ibin++){
    resGr=Test.GetResErr(nobs+ibin);	//bins final results and errors
    OutHisto<<setw(wd)<< ibin <<setw(wd)<< resGr[0] <<setw(wd)<< resGr[1]<<endl;
}
OutHisto.close();

return 0;
}
