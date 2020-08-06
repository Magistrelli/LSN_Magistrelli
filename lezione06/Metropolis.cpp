#include "Metropolis.h"

Metropolis::Metropolis(): Experiment(),UseGauss(0) {}

Metropolis::Metropolis(Random* rnd): Experiment(rnd),HNow(0),HNew(0),Beta(0) {
    UseGauss=0;			//default: uniform T
    Ending=0;
}

//effective single step of Markov Chain, std move with uniform or Gaussian T
//(ip not important if use StdXnew)
void Metropolis::StdMove(DataVett Xnew, int ip){
    bool ifMove=StdIfMove(Xnew,ip);
    if(ifMove==1){			//new point accepted
	if(ip==-1) {*X=Xnew;}
	else {				//Xnew has only DimSp component
	    for(int k=0;k<DimSp;++k) {X->SetComp(DimSp*ip+k,Xnew.GetComp(k));}
	}
    }
}
//check if move accepted or not
bool Metropolis::StdIfMove(DataVett Xnew, int ip){
    double q,A=1.,sort;
    bool ifMove=0;
    q=qRatio(Xnew,ip);			//T(x|y)=T(y|x) => q=p(Xnew)/p(Xold)
    if(q<1.)	{A=q;}
    sort=Rnd->Rannyu();
    if(sort<=A){			//new point accepted
	ifMove=1;
	accepted++;
    }
    tries++;
    return ifMove;
}

//new configuration's standard generation
DataVett Metropolis::StdXnew(double param) const{
    DataVett Xnew(TotD);
    for(int i=0; i<nPart; ++i){
      for(int j=0; j<DimSp; j++){
	if(UseGauss==1){//Xnew sorted with Gaussian T, mu=0.,sigma=param
	    Xnew.DefComp(DimSp*i+j,X->GetComp(DimSp*i+j)+Rnd->Gauss(0.,param));
	} else {	//Xnew sorted with uniform T, width param
	    Xnew.DefComp(DimSp*i+j,X->GetComp(DimSp*i+j)+Rnd->Rannyu(-param,param));
	}
      }
    }
    return Xnew;
}


//simulated annealing algorithm to find minimum of AnnEnergy
//args describe annealing schedule (beta_i,n_i)
void Metropolis::SimAnnealing(const DataVett beta, const DataVett bStep){
    if(beta.GetDim()!=bStep.GetDim()){
	cerr<<endl<<"Attention! Annealing schedule's dim not compatible!";
	exit(-1);
    }
    ofstream Null[1];
    HNow=AnnEnergy();		//HNow initialized here

    cout << "Simulated Annealing method:" << endl;
    for(int i=0; i<int(beta.GetDim()); ++i){
	Beta=beta.GetComp(i);	//actual temperature
	cout << setprecision(2) << "  sampling T ~ " << 1./Beta;
	Reset(1);		//reset accumulators
	for(int istep=0; istep<bStep.GetComp(i); ++istep) {Move();}
	PrintAcc();
    }//now we are in the minimum region, we can do final blocking method meas:
    Ending=1;
    cout << endl << "Doing final SA exstimation:" << endl;
    DoSampling(Null,Null[0],0);
    cout << "  Acceptance: " << Acceptance();   
}

//qRatio for Simulated Annealing
double Metropolis::AnnqRatio() const {
    if(HNew>HNow)	{return exp(-Beta*(HNew-HNow));}//accepting prob
    else		{return 1.;}			//move accepted
}
