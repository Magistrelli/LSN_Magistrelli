#include "MolDyn_NVE.h"

//constructors
MolDyn::MolDyn(): Experiment(),Temp(0),Rho(0),Vol(0),Box(0),Rcut(0),Delta(0),BinSize(0)	{Xold=NULL,V=NULL;}

MolDyn::MolDyn(int nblk, Random* rnd): Experiment(rnd) {//Prepare all stuff for the simulation
    double sumv[3]={0.,0.,0.},vel=0.,in;
    ifstream ReadInput,ReadConf,ReadOld;

    DimSp=3;
    m_Props=n_Props;
    ReadInput.open("input.dat");	//read input
    ReadInput >> Temp;
    ReadInput >> nPart;
    TotD=DimSp*nPart;
    ReadInput >> Rho;
    Vol=(double)nPart/Rho;
    Box=pow(Vol,1./3.);
    BinSize=(Box/2.)/double(nBins);
    ReadInput >> Rcut;
    ReadInput >> Delta;
    nBlk=nblk;
    ReadInput >> nStep;
    if (nStep%nBlk!=0){
	cerr << endl << "==================================" << endl;
	cerr <<" Attention! (nStep)%(nBlk)!=0 !!! " << endl;
	cerr <<"==================================" << endl << endl;
    }
    nStep/=nBlk;
    ReadInput >> iPrint;
    ReadInput >> _restart;	//0 for false, 1 for true
    ReadInput.close();
    PrintInfo();

    //Read starting configuration
    InitVett();			//init X, Walker, ...
    ReadConf.open("config.0");
    cout << "Read initial configuration from file config.0 " << endl;
    for (int i=0; i<TotD; ++i){
	ReadConf >> in;
	X->DefComp(i,Pbc(in*Box));
    }
    ReadConf.close();

    Xold=new DataVett(TotD);	//(1x,1y,1z,2x,2y,2z,...,nPart_z)
    V=new DataVett(TotD);
    //ReadOld if restart, velocity exstimation from actual and new configuration, new starting point
    if(_restart==1){
      ReadOld.open("old.0");
      if (ReadOld.fail()) {
        cout << endl << "File old.0 doesn't exist!" << endl << "\tSimulation will compute old position via starting velocity exstimation" << endl;
        _restart=0;
      }
      else{
	cout << "Read old configuration from file old.0 " << endl << endl;
	for (int i=0; i<TotD; ++i){
	    ReadOld >> in;
	    Xold->DefComp(i,Pbc(in*Box));
	}
	Restart();
      }
      ReadOld.close();
    }
    //if not restart, compute old via starting velocity exstimation
    if(_restart==0){
	cout << "Prepare random starting velocities with center of mass velocity equal to zero " << endl << endl;
	for (int i=0; i<nPart; ++i){
	    for (int j=0; j<DimSp; j++){
		V->DefComp(DimSp*i+j,Rnd->Rannyu(-0.5,0.5));//random vel
		sumv[j]+=V->GetComp(DimSp*i+j);//total velocity, each component
	    }
	}
	for (int j=0; j<DimSp; ++j)	{sumv[j]/=(double)nPart;}
	for (int i=0; i<nPart; ++i){	//each component of total V at 0
	    for (int j=0; j<DimSp; j++)	{V->SumComp(DimSp*i+j,-sumv[j]);}
	}
	for(int i=0; i<TotD; ++i) {vel+=pow(V->GetComp(i),2);}
	(*V)*=sqrt(3*Temp/(vel/double(nPart)));	//velocity rescaled with T
	*Xold=Pbc(*X-Delta*(*V));
    }
}

//distructor
MolDyn::~MolDyn(){
    delete Xold;
    delete V;
}


//restart simulation for warm-up
void MolDyn::Restart(){
    double Tnow,vel=0.;
    V->SetUsed(0);		//re-initialize V
    Move(X,Delta);//one step of Verlet algorithm but compute V(t+dt/2)
    for(int i=0; i<TotD; ++i) {vel+=pow(V->GetComp(i),2);}
    Tnow=vel/(3.*double(nPart));//exstimate of T(t+dt/2)
    (*V)*=sqrt(Temp/Tnow);	//velocity rescaled with temperature
    *Xold=Pbc(*X-Delta*(*V));	//new x(t) to match temp
    Reset(1);			//accumulated variables to 0
}


//Compute forces as -Grad_ip V(r)
double MolDyn::Force(const int ip, const int idir) const{
    double f=0.,dvec[3],dr,a,b;
    for (int i=0; i<nPart; ++i){
	if(i!=ip){
	    dr=0.;
	    for (int j=0; j<DimSp; ++j){
		a=X->GetComp(DimSp*ip+j),b=X->GetComp(DimSp*i+j);
		dvec[j]=Pbc(a-b);		// distance ip-i in pbc
		dr+=dvec[j]*dvec[j];
	    }
	    dr=sqrt(dr);
	    if(dr<Rcut) {f+=dvec[idir]*(2./pow(dr,14)-1./pow(dr,8));} // -Grad_ip V(r)
	}
    }
    return f*24.;
}


//function used in InstantValues
void MolDyn::WriteInstant(int istep, ofstream* OutRes) const{
    const int wd=12;
    for(int i=0; i<iGr; i++) {OutRes[i]<<setw(wd)<<istep<<setw(wd)<<Walker->GetComp(i)/double(nPart)<<endl;}
}


//Move particles with Verlet algorithm (compute V(t+jump/2) with new and Y positions
//with this method I define stdMove in MolDyn.h
void MolDyn::Move(DataVett* Y, double jump){
    DataVett Xnew;
    Xnew=2.*(*X)-*Xold;
    for(int i=0; i<nPart; ++i){		//Verlet integration scheme
	for (int j=0; j<DimSp; ++j) {Xnew.SumComp(DimSp*i+j,Force(i,j)*pow(Delta,2));}
    }
    Xnew=Pbc(Xnew);
    *V=Pbc(Xnew-*Y)/jump;
    *Xold=*X;
    *X=Xnew;
}

//Properties measurement
void MolDyn::Measure(){
    double dx,dr,bin_m,v=0.,nowKin=0.,nowPot;
    for (int i=0; i<nPart-1; ++i){	//cycle over pairs of particles
      for (int j=i+1; j<nPart; ++j){
	dr=0.;
	for (int k=0; k<DimSp; ++k){	//use old configurations [old = r(t)] to be compatible with EKin which uses v(t) => EPot should be computed with r(t)
	    dx=Pbc(Xold->GetComp(DimSp*i+k)-Xold->GetComp(DimSp*j+k));
	    dr+=dx*dx;
	}
	dr = sqrt(dr);

	for(int ibin=nBins-1; ibin>=0; --ibin){//update of the histogram of g(r)
	    bin_m=ibin*BinSize;		//bin's left extreme
	    if((dr>bin_m)&&(dr<bin_m+BinSize)){
	      Walker->SumComp(iGr+ibin,2);
	      ibin=-1;			//exit if bin finded
	    }
	}
	if(dr<Rcut)	{v+=1./pow(dr,12)-1./pow(dr,6);}//Potential energy
      }
    }
    for(int i=0; i<TotD; ++i) {nowKin+=pow(V->GetComp(i),2);}
    nowPot=4.*v,nowKin*=0.5;
    Walker->SumComp(iV,nowPot);
    Walker->SumComp(iK,nowKin);
    Walker->SumComp(iE,nowPot+nowKin);
    Walker->SumComp(iT,2./3.*nowKin);
}
//Print results for current block
void MolDyn::Averages(int iblk, ofstream *OutRes){
    double gAv,bin_m,bin_M,dVr;
    DataVett stima(n_Props);
    const int wd=12;
    
    for(int i=0; i<iGr; ++i)	{stima.DefComp(i,Walker->GetComp(i)/double(nStep*nPart));}
    for(int ibin=0; ibin<nBins; ++ibin){
        gAv=Walker->GetComp(iGr+ibin)/double(nStep);
	bin_m=ibin*BinSize;		//bin's left extreme
	bin_M=bin_m+BinSize;		//bin's right extreme
	dVr=4./3.*M_PI*(pow(bin_M,3)-pow(bin_m,3));
	stima.DefComp(iGr+ibin,gAv/(Rho*nPart*dVr));//g(r)
    }

    DoAverages(stima,iblk);	//do blocks progressive averages
    m_Props=iGr;		//write progressive res and err only for V,K,E,T
    WriteAverages(stima,iblk,OutRes);//write output
    m_Props=n_Props;		//(modified m_Props only for WriteAverages)

    if(iblk==1)	{OutRes[iGr] << setw(wd) << "block #" << setw(2*wd) << "bins" << endl << endl;}
    OutRes[iGr]<< setw(wd) << iblk << setw(2*wd);
    for(int i=0; i<nBins; i++)	{OutRes[iGr]<< stima.GetComp(iGr+i) << setw(wd);}
    OutRes[iGr] << endl;	//here I wrote only single blocks results
}


//Pbc applied to an entire DataVett
DataVett MolDyn::Pbc(const DataVett& Y) const {
    const int dim=Y.GetDim();
    DataVett r(dim);
    for (int i=0; i<dim; ++i)	{r.DefComp(i,Pbc(Y.GetComp(i)));}
    return r;
}


//Write final configuration, new and old positions
void MolDyn::ConfFinalPlus() const{
    ofstream WriteOld;
    cout << "Print final old configuration to file old.final " << endl;
    WriteOld.open("old.final");
    for (int i=0; i<nPart; ++i){
	for (int j=0; j<DimSp; j++) {WriteOld << Xold->GetComp(DimSp*i+j)/Box << "   ";}
	WriteOld << endl;
    }
    WriteOld.close();
    ConfFinal();
}


//Print simulation's information and parameters
void MolDyn::PrintInfo() const {
    cout << endl << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  "<< endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;
    
    cout << "Number of particles = " << nPart << endl;
    cout << "Density of particles = " << Rho << endl;
    cout << "Volume of the simulation box = " << Vol << endl;
    cout << "Edge of the simulation box = " << Box << endl;

    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << Delta << endl;
    cout << "Number of total steps = " << nStep*nBlk << endl << endl;
}
