#include "Ising.h"

int main (int argc, char *argv[]){
if (argc!=3) {
    cerr << "Usage: " << argv[0] << " <Input Temperature parameter> <restarting index (#restart-1)>" << endl;
    return -1;
}

//Variables
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");//initialization random_gen from files
Ising Test(rnd);
int nT=20,nObs=4;	//number of different temperatures -1, obs sampled
double T_in=(2.-0.5)/double(nT)*atof(argv[1])+0.5;//input temperature
double *reserr=new double[2];//single T_in result and error couple
ofstream OutRes[nObs],OutBurn,OutFinal,Null[1];
string *namefolder,*nameFile,*metro;
if(Test.GetMetro()==1)	{metro=new string("metro");}
else			{metro=new string("gibbs");}
const int wd=20;

//point 1.
namefolder=new string("results");
if(Test.Geth()!=0) {*namefolder+="_h";}

Test.SetTemp(T_in);
Test.PrintInfo();
if(atoi(argv[2])==0){		//first simulation's run
    Test.SetBurn(int(5e3));
    OutBurn.open(*namefolder+"/burnIn."+*metro+"_"+to_string(atoi(argv[1]))+".out");
    Test.BurnIn(OutBurn,1);
    OutBurn.close();
}
Test.DoSampling(Null,Null[0],0);

nameFile=new string(*namefolder+"/config."+*metro+"."+to_string(atoi(argv[1]))+".final");
cout << "Print final configuration to file " << *nameFile << endl << endl;
OutFinal.open(*nameFile);
Test.ConfFinal(OutFinal);	//save final config
OutFinal.close();

OutRes[0].open(*namefolder+"/output."+*metro+".ene."+to_string(atoi(argv[2])),ios::app);
OutRes[1].open(*namefolder+"/output."+*metro+".heat."+to_string(atoi(argv[2])),ios::app);
OutRes[2].open(*namefolder+"/output."+*metro+".mag."+to_string(atoi(argv[2])),ios::app);
OutRes[3].open(*namefolder+"/output."+*metro+".chi."+to_string(atoi(argv[2])),ios::app);
for(int iObs=0; iObs<nObs; iObs++) {
    reserr=Test.GetResErr(iObs);
    OutRes[iObs] << setw(wd) << T_in << setw(wd) << reserr[0] << setw(wd) << reserr[1] << endl;
}
for(int i=0; i<nObs; i++){OutRes[i].close();}

cout << "Temperature: " << T_in << endl;
cout << "Mean acceptance rate: " << Test.Acceptance() << endl << endl;
cout << "----------------------------" << endl << endl;

return 0;
}
