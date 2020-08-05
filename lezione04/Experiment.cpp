#include "Experiment.h"

//EXPERIMENT
Experiment::Experiment(): nPart(0),DimSp(0),TotD(0),nBlk(0),nStep(0),nBurn(0),accepted(0), tries(0),TotA(0),iPrint(0),_restart(0),m_Props(0){
    Rnd=NULL;
    X= NULL;
    Walker=NULL;
    GlobAv=NULL,GlobAv2=NULL;
    Err= NULL;
}

Experiment::Experiment(Random* rnd) {
    nPart=0;
    DimSp=1;				//default: 1-D space
    nBlk=0,nStep=0;
    nBurn=0;
    accepted=0,tries=0;			//(p.e. Metropolis' 50% criterion)
    TotA=0.;
    Rnd=rnd;
    _restart=0;				//deafaul: false
    X=NULL;
    m_Props=0;
    Walker=NULL;
    GlobAv=NULL,GlobAv2=NULL;
    Err= NULL;
}

Experiment::~Experiment(){ //don't delete Rnd, maybe useful also for main
    delete X;
    delete Walker;
    delete GlobAv;
    delete GlobAv2;
    delete Err;
}

//DataVett initialization
void Experiment::InitVett(){
    X=new DataVett(TotD);	//(1x,1y,1z,2x,2y,2z,...,nPart_z)
    Walker=new DataVett(m_Props);
    GlobAv=new DataVett(m_Props);
    GlobAv2=new DataVett(m_Props);
    Err=new DataVett(m_Props);
    Reset(1);
}


//Reset block averages
void Experiment::Reset(int iblk){
    DataVett Null(m_Props);
    Null.SetUsed(m_Props);
    if(iblk==1){
	TotA=0.;
	*GlobAv=Null;
	*GlobAv2=Null;
	*Err=Null;
    }
    *Walker=Null;
    tries=0;
    accepted=0;
}

//Re-initialize
void Experiment::Restart(const DataVett start){
    for(int i=0; i<TotD; ++i){X->SetComp(i,start.GetComp(i));}
    Reset(1);
}


//if progressive evolution of res and err not needed, can extract final res
double* Experiment::GetResErr(int iObs) const {
    double *out = new double[2];
    out[0]=GlobAv->GetComp(iObs)/double(nBlk);
    out[1]=Err->GetComp(iObs);
    return out;
}


//Effective computation observables' averages
void Experiment::DoAverages(DataVett stima, int iblk) const {
    *GlobAv+=stima;		//cumulative sum of blocks results
    *GlobAv2+=stima*stima;	//sum for variance
    for(int i=0; i<m_Props; ++i) {Err->SetComp(i,Error(GlobAv->GetComp(i),GlobAv2->GetComp(i),iblk));}
}
//Progressive result-error files printing
void Experiment::WriteAverages(DataVett stima, int iblk, ofstream* OutRes) const {
    const int wd=20;
    for(int i=0; i<m_Props; i++) {OutRes[i]<<setprecision(8)<<setw(wd)<<iblk<<setw(wd)<<stima.GetComp(i)<<setw(wd)<<GlobAv->GetComp(i)/(double)iblk<<setw(wd)<<Err->GetComp(i)<<endl;}
}


//print last nnblk blocks' configurations in .config
void Experiment::WriteConf(ofstream& OutConf, const int iblk, const int istep, int nnblk) const {
    if (!OutConf.fail()){
      if(iblk>(nBlk-nnblk)){	//Only for the last nnblk blocks
	for(int i=0; i<nPart; ++i){
	    for(int j=0; j<DimSp; ++j)	{OutConf << X->GetComp(DimSp*i+j) << " ";}
	    OutConf << endl;
	}
	if(istep==nStep){cout<<endl<<"Print configurations to file .config";}
      }
    }
}

//Write configuration in .xyz format
void Experiment::ConfXYZ(int nconf) const{
    ofstream WriteXYZ;
    WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
    WriteXYZ << nPart << endl;
    WriteXYZ << "This is only a comment!" << endl;
    for(int i=0; i<nPart; ++i){
	WriteXYZ << "LJ  ";
	for(int j=0; j<DimSp; ++j) {WriteXYZ << X->GetComp(DimSp*i+j) <<" ";}
	WriteXYZ << endl;
    }
    WriteXYZ.close();
}

//Write final config (deafaul file)
void Experiment::ConfFinal() const{
    ofstream Out;
    cout<<"Print final configuration to file config.final "<<endl<<endl;
    Out.open("config.final");
    ConfFinal(Out);
    Out.close();
}
//Write final configuration
void Experiment::ConfFinal(ofstream& Out) const{
    for(int i=0; i<nPart; ++i){
	for(int j=0; j<DimSp; ++j) {Out << X->GetComp(DimSp*i+j) << "   ";}
	Out << endl;
    }
    Rnd->SaveSeed();
}


//Core of the experiment simulation
//if conf==1 write last block's configurations and print progresses
void Experiment::DoSampling(ofstream* OutRes, ofstream& OutConf, int nnconf){
    for(int iblk=1; iblk<=nBlk; ++iblk){	//Simulation
	if ((nBlk>5)&(nBlk<=20)){cout <<"doing block "<< iblk << "/" << nBlk <<endl;}
	else if (iblk%5==0)	{cout <<"doing block "<< iblk << "/" << nBlk <<endl;}
	Reset(iblk);				//Reset block averages
	for(int istep=1; istep<=nStep; ++istep){
	    Move();
	    Measure();
	    if (nnconf>0){ //write actual conf if in last nnconf blocks
		WriteConf(OutConf,iblk,istep,nnconf);
		OutConf << endl << endl;
	    }
	}
	Averages(iblk,OutRes);		//Print results for current block
	TotA+=double(accepted)/double(tries);
    }
}


//Move particles for nBurn steps and write instant measure and single configurations for the first nobs obs
void Experiment::InstantValues(ofstream* OutRes, int nobs, int iout, int iconf) {
    int nconf=1;
    DataVett Null(m_Props);
    Null.SetUsed(m_Props);

    for(int i=0; i<nobs; i++) {OutRes[i] << nBurn << endl << endl;}
    for(int istep=1; istep<=nBurn; ++istep){
	Move();
	if(istep%iout==0){
	    *Walker=Null;	//reset every time, instant results
	    Measure();		//Properties measurement
	    WriteInstant(istep,OutRes);
	}
	if (iconf>0) {	if (istep%iconf==0){ //if iconf<=0 doesn't write conf
			    ConfXYZ(nconf);
			    nconf++;
			}
	}
	if(istep%iPrint==0) {cout<<"Number of time-steps: "<< istep << endl;}
    }
    if ((accepted!=0)&(tries!=0)) {cout << endl << "Mean Acceptance Rate = " << double(accepted)/double(tries) << endl;} //print only if used
    cout << endl;
}


//burn-in (write in default file)
void Experiment::BurnIn(){
    ofstream OutBurn;
    cout << "Print burn-in to file burn_in.out " << endl << endl;
    OutBurn.open("burn_in.out");
    BurnIn(OutBurn,1);
    OutBurn.close();
}
//burn-in and write istant values of Walker[0] to visualize burn-in
void Experiment::BurnIn(ofstream& OutBurn, bool out){
    const int iprint=1, wd=12;
    if (out==1)	{OutBurn << nBurn << endl << endl;}
    for(int k=1; k<=nBurn; k++){
	Move();
	if ((out==1)&(k%iprint==0)) {
	    Walker->SetComp(0,0.);	//reset every time, instant results
	    Measure();
	    OutBurn<<setw(wd)<< k <<setw(wd)<< Walker->GetComp(0) <<endl;
	}
    }
}





//EXTERNAL FUNCTIONS

//statistical uncertanties for Blocking Method
double Error(double sum, double sum2, int iblk){
    double ave,av2,sigma;
    ave=sum/double(iblk);
    av2=sum2/double(iblk);
    sigma=sqrt((av2-ave*ave)/(iblk-1.));
    if(iblk==1)	{return 0.;}
    else	{return sigma;}
}
