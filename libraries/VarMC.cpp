#include "VarMC.h"

//PSITRIAL
#define MuS X	//for 2-gaussian PsiTrial X = parameters (mu,sigma)

//constructor
PsiTrial::PsiTrial(Random* rnd): Metropolis(rnd),Delta(0){
    nPart=1,DimSp=2;		//1D function, 2D parameters' space
    TotD=DimSp*nPart;
    m_Props=n_Props;
    Metro=NULL;			//will be initialized only for SimAnnealing
    InitVett();
    MuS->DefComp(0,1.);		//start in mu = 1
    MuS->DefComp(1,1.);		//start in sigma = 1
}
PsiTrial::PsiTrial(Random* rnd, double delta): PsiTrial(rnd) {Delta=delta;}

//evaluate PsiTrial (here superposition of two gaussians, but generalizable)
//PsiT modification requires also Der2 modification
double PsiTrial::Eval(double x) const {
    double mu=MuS->GetComp(0),sig2=2.*pow(MuS->GetComp(1),2);
    return exp(-pow(x-mu,2)/sig2)+exp(-pow(x+mu,2)/sig2);
}
//second derivative
double PsiTrial::Der2(double x) const {
    double mu=MuS->GetComp(0),sig2=pow(MuS->GetComp(1),2);
    double xpm2=pow(x+mu,2),xmm2=pow(x-mu,2);
    return exp(-xpm2/(2.*sig2))*(exp(2.*mu*x/sig2)*(xmm2-sig2)+xpm2-sig2)/pow(sig2,2);
}


//check and save new config
void PsiTrial::Move() {
    DataVett parOld=*MuS;	//save actual parameters
    *MuS=StdXnew(Delta);	//set new trial parameters for Metro sampling
    HNew=AnnEnergy();		//measure <H>new, no output    
    if(StdIfMove(DataVett(0),-1)){//move accepted;
	HNow=HNew;		//save new stored energy (don't resample)
    } else {			//move rejected
	*MuS=parOld;		//restore old parameters
    }
}

//Evaluate local energy and Psi(r) of the actual configuration
void PsiTrial::Measure(){
    Walker->SumComp(iM,MuS->GetComp(0));
    Walker->SumComp(iS,MuS->GetComp(1));
}
//Print results for current block
void PsiTrial::Averages(int iblk, ofstream *OutRes){
    DataVett stima(n_Props);
    stima.DefComp(iM,Walker->GetComp(iM)/double(nStep));
    stima.DefComp(iS,Walker->GetComp(iS)/double(nStep));
    DoAverages(stima,iblk);
}


//energy eval for Simulated Annealing; here is Metro's <H> exstimation
double PsiTrial::AnnEnergy() const {
    ofstream Null[1];
    Metro->BurnIn(Null[0],0);		//rapid (start from 0) burn-in
    Metro->DoSampling(Null,Null[0],0);	//measure initial <H>, no output
    return (Metro->GetResErr(0))[0];	//HNow initialized here
}
//print psi and last VarMC acceptance parameters for every temperature
void PsiTrial::PrintAcc() const {
    cout << "\tAcceptance:\tpsi = "<< setprecision(6) << double(accepted)/double(tries)<<"\t\tMetro (last) = "<< Metro->Acceptance() <<endl;
}


#undef MuS




//VARMC

//constructor
VarMC::VarMC(PsiTrial* psi, double xl, double del, Random* rnd): Metropolis(rnd){
    nPart=1,DimSp=1;	//1 particle, 1D space
    TotD=DimSp*nPart;
    m_Props=n_Props;
    iPrint=int(1e7);
    PsiT=psi;
    Delta=del;
    Xlim=xl;
    BinSize=2.*Xlim/double(nBins);

    InitVett();		//X = PsiTrial's parameters (mu,sigma)
    X->DefComp(0,0.);	//set x0=4. to see burn-in in the first 10 points
}


//Hamiltonian applied to PsiTrial (H=T+V, T=p2/2m) divided for PsiTrial
double VarMC::Eloc(double x) const
    {return -(hbar*hbar)*PsiT->Der2(x)/(2.*Mass)/PsiT->Eval(x)+VExt(x);}


//Evaluate local energy and Psi(r) of the actual configuration
void VarMC::Measure(){
    double bin_m,x=X->GetComp(0);
    for(int ibin=0; ibin<nBins; ++ibin){//update of the histogram of psi(r)
	bin_m=-Xlim+ibin*BinSize;	//bin's left extreme
	if((x>bin_m)&(x<bin_m+BinSize)){
	    Walker->SumComp(iF+ibin,1);
	    ibin=nBins;			//bin found
	}
    }
    Walker->SumComp(iE,Eloc(x));
}
//Print results for current block
void VarMC::Averages(int iblk, ofstream *OutRes){
    double fAv[nBins],norm=0.;
    DataVett stima(n_Props);
    
    stima.DefComp(iE,Walker->GetComp(iE)/double(nStep*nPart));
    for(int ibin=0; ibin<nBins; ++ibin){	//psi(x)
        fAv[ibin]=Walker->GetComp(iF+ibin)/double(nStep*nPart);
	norm+=fAv[ibin];			//normalization
    }
    norm*=BinSize;
    for(int ibin=0;ibin<nBins;++ibin) {stima.DefComp(iF+ibin,fAv[ibin]/norm);}
    
    DoAverages(stima,iblk);
    m_Props=iF;			//write progressive res and err only for <H>
    WriteAverages(stima,iblk,OutRes);
    m_Props=n_Props;		//(modified m_Props only for WriteAverages)
}

//T(x|y)=T(y|x) => q=p(Xnew)/p(Xold)
double VarMC::qRatio(const DataVett Xnew, int ip) const{
    if (Xnew.GetDim()!=1){
	cerr << "Attention! 1D single-particle problem! Given Xnew of a " << Xnew.GetDim() << "D space" << endl;
	exit(-1);
    }
    return pow(PsiT->Eval(Xnew.GetComp(0))/PsiT->Eval(X->GetComp(0)),2);
}
