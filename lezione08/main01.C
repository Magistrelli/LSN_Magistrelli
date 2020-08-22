#include "VarMC.h"

int main (int argc, char *argv[]){
if (argc!=3) {
    cerr << "Usage: " << argv[0] << " <mu (PsiT par)> <sigma (PsiT par)>" << endl;
    return -1;
}

//Variables
DataVett MuS(2);			//chosen mu and sigma
double xlim=5.,delta=2.5;		//studied x in [-5,5], unif-T width
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");	//initialization random_gen from files
PsiTrial* psi=new PsiTrial(rnd);	//trial wave functions Psi(mu,sigma)
VarMC Test(psi,xlim,delta,rnd);
double* res=new double[2];
ofstream OutBurn,OutRes[1],Null;	//obs is energy <H>

//simulation
cout<<endl;
for(int i=0; i<2; ++i)	{MuS.DefComp(i,atof(argv[i+1]));}
psi->SetPar(MuS);
Test.SetBlkSt(100,int(1e6));

Test.SetBurn(int(1e3));
OutBurn.open("results/trial."+string(argv[1])+"_"+string(argv[2])+".ene.burn");
Test.BurnIn(OutBurn,1);
OutBurn.close();

OutRes[0].open("results/trial."+string(argv[1])+"_"+string(argv[2])+".ene.out");
Test.DoSampling(OutRes,Null,0);
OutRes[0].close();

res=Test.GetResErr(0);
cout<<endl <<"Mean Acceptance Rate = " << Test.Acceptance();
cout<<endl <<"Ground state energy: "<< res[0] <<" +- "<< res[1] <<endl<<endl;

return 0;
}
