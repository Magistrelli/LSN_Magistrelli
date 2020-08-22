#include "VarMC.h"

int main (int argc, char *argv[]){

//Variables
double xlim=5.,delVar=2.6,delPsi=0.15;	//studied x in [-5,5], unif-Ts width
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");	//initialization random_gen from files
PsiTrial* psi=new PsiTrial(rnd,delPsi);	//trial wave functions Psi(mu,sigma)
VarMC* Test=new VarMC(psi,xlim,delVar,rnd);
const int nbeta=20,nbstep=200;		//#annealing_beta, #step for each beta
const int psiBlk=10,psiStp=int(1e3);	//#blk, #stp for psi params final meas
const int Hblk=100,Hstp=int(1e6);	//#blk, #stp for <H> final measure
const double bmax=1./0.01,bmin=1./10.;	//min and max annealing beta
DataVett beta(nbeta+1),bStep(nbeta+1);	//annealing schedule
double* mus=new double[2];		//res from annealing
double* res=new double[2];		//res for <H>
DataVett MuS(2),MuSerr(2);		//optimized parameters and errors
const int nobs=1,nbins=Test->GetNBins();
ofstream OutRes[nobs],OutHisto,Output,Null;//obs is energy <H>
const int wd=20;

//simulation
cout<<endl;
Test->SetBurn(100);		//rapid (start from 0) burn-in
Test->SetBlkSt(1,int(5e3));	//no uncentainty here for <H> exstimation
psi->SetBlkSt(psiBlk,psiStp);
for(int i=0; i<=nbeta; ++i){	//define annealing schedule
    beta.DefComp(i,bmin+i*(bmax-bmin)/nbeta);
    bStep.DefComp(i,nbstep);
}
psi->SetMetro(Test);
psi->SimAnnealing(beta,bStep);

for(int i=0; i<2; ++i){
    mus=psi->GetResErr(i);
    MuS.DefComp(i,mus[0]);	//psi parameters values
    MuSerr.DefComp(i,mus[1]);	//and statistical uncertainty
}
delete mus;

cout<<endl<<endl<< "Exstimation of <H> with optimized mu and sigma:" << endl;
psi->SetPar(MuS);		//set optimized parameters in psi
Test->SetBlkSt(Hblk,Hstp);
Test->BurnIn(Null,0);
OutRes[0].open("results/var.ene.out");
Test->DoSampling(OutRes,Null,0);//compute final estimation of <H>
OutRes[0].close();
res=Test->GetResErr(0);
cout << "  Mean Acceptance Rate = " << Test->Acceptance() << endl << endl;
cout<<"Print energy results to file results/var.ene.out"<<endl;

Output.open("results/psiTrial.optimized.out");
Output << MuS.GetComp(0) << endl << MuS.GetComp(1) << endl << res[0] << endl << MuSerr.GetComp(0) << endl << MuSerr.GetComp(1) << endl << res[1] << endl << endl << "Minimum parameters:\n  mu;\n  sigma;\n  <H>;\n  mu_uncertainty;\n  sigma_uncertainty;\n  <H> uncertainty;";
Output.close();

OutHisto.open("results/var.psi.out");
OutHisto << Hstp << setw(wd) << Test->GetBinWidth() << endl << endl;
for(int ibin=0; ibin<nbins; ibin++){
    res=Test->GetResErr(nobs+ibin);
    OutHisto<<setw(wd)<< ibin <<setw(wd)<< res[0] <<setw(wd)<< res[1] << endl;
}
OutHisto.close();
cout<<"Print psi(x) for histogram to file results/var.psi.out"<<endl<<endl;

delete psi;
delete Test;
delete res;
return 0;
}
