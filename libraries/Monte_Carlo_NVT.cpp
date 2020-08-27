#include "Monte_Carlo_NVT.h"

//constructors
CanonicEns::CanonicEns(Random* rnd): Metropolis(rnd){
    ifstream ReadInput,ReadConf;
    double in;

    DimSp=3;
    m_Props=n_Props;
    iPrint=int(2.5e4);
    ReadInput.open("input.dat");
    ReadInput >> Temp;
    Beta = 1./Temp;
    ReadInput >> nPart;
    TotD=DimSp*nPart;
    ReadInput >> Rho;
    Vol = (double)nPart/Rho;
    Box = pow(Vol,1./3.);
    BinSize = (Box/2.)/double(nBins);
    ReadInput >> Rcut;
    ReadInput >> Delta;
    Delta/=2.;
    ReadInput >> nBlk;
    ReadInput >> nStep;
    ReadInput.close();
    //Tail corrections for potential energy and pressure
    Vtail=(8.*M_PI*Rho)/(9.*pow(Rcut,9))-(8.*M_PI*Rho)/(3.*pow(Rcut,3));
    Ptail=(32.*M_PI*Rho)/(9.*pow(Rcut,9))-(16.*M_PI*Rho)/(3.*pow(Rcut,3));
    PrintInfo();

    //Read initial configuration
    InitVett();
    ReadConf.open("config.0");
    cout<<endl<<"Read initial configuration from file config.0"<<endl<<endl;
    for(int i=0; i<TotD; ++i){
	    ReadConf >> in;
	    X->DefComp(i,Pbc(in*Box));
	}
    ReadConf.close();

    Measure();	//Evaluate potential energy and virial of the initial config
    cout << "Initial values (with tail corrections):" << endl;
    cout << "Potential energy = "<<Epot(Walker->GetComp(iV)) << endl;
    cout << "Virial           = "<<Virial(Walker->GetComp(iW)) << endl;
    cout << "Pressure         = "<<Pressure(Walker->GetComp(iW))<<endl<<endl;
}


void CanonicEns::WriteInstant(int istep, ofstream* OutRes) const {
    const int wd=20;
    OutRes[iV]<<setw(wd)<<istep<<setw(wd)<< Epot(Walker->GetComp(iV)) <<endl;
    OutRes[iW]<<setw(wd)<<istep<<setw(wd)<< Pressure(Walker->GetComp(iW)) <<endl;
}


//prob ratio, proposed new config by changing only the ip particle position
double CanonicEns::qRatio(const DataVett Xnew, int ip) const{
    double dx,dx_new,dr,dr_new,ene=0.,ene_new=0.;
    for (int i=0; i<nPart; ++i){
      if(i!=ip){
        dr=0.,dr_new=0.;
	for(int j=0; j<DimSp; ++j){	// distance ip-i in pbc
	    dx=Pbc(X->GetComp(DimSp*ip+j)-X->GetComp(DimSp*i+j));
	    dx_new=Pbc(Xnew.GetComp(j)-X->GetComp(DimSp*i+j));
	    dr+=dx*dx,dr_new+=dx_new*dx_new;	    
	}
	dr=sqrt(dr),dr_new=sqrt(dr_new);
	if(dr<Rcut)	{ene+=1./pow(dr,12)-1./pow(dr,6);}
	if(dr_new<Rcut)	{ene_new+=1./pow(dr_new,12)-1./pow(dr_new,6);}
      }
    }
    ene*=4.;
    ene_new*=4.;
    return exp(Beta*(ene-ene_new));
}


//one step of Markov Chain, try to change nPart random particles
void CanonicEns::Move(){
    int o;
    DataVett Xnew(DimSp);
    Xnew.SetUsed(DimSp);
    for(int i=0; i<nPart; ++i){		//108 tries every MC step
	o=int(Rnd->Rannyu()*nPart);	//Random particle (0<=o<=npart-1)
	for(int j=0; j<DimSp; ++j)	{Xnew.SetComp(j,Pbc(X->GetComp(DimSp*o+j)+Rnd->Rannyu(-Delta,Delta)));}
	StdMove(Xnew,o);
    }
}

//Evaluate potential energy, virial and g(r) of the actual configuration
void CanonicEns::Measure(){
    double dx,dr,bin_m,v=0.,w=0.;
    for (int i=0; i<nPart-1; ++i){	//cycle over pairs of particles
      for (int j=i+1; j<nPart; ++j){	//j>i
	dr=0.;
	for(int k=0; k<DimSp; ++k){	//distance i-j in pbc
	    dx=Pbc(X->GetComp(DimSp*i+k)-X->GetComp(DimSp*j+k));
	    dr+=dx*dx;
	}
	dr=sqrt(dr);

	for(int ibin=nBins-1; ibin>=0; --ibin){//update of the histogram of g(r)
	    bin_m=ibin*BinSize;		//bin's left extreme
	    if((dr>bin_m)&&(dr<bin_m+BinSize)){
	      Walker->SumComp(iGr+ibin,2);
	      ibin=-1;			//exit if bin finded
	    }
	}
	if(dr<Rcut){		//contribution to energy and virial
	    v+=1./pow(dr,12)-1./pow(dr,6);
	    w+=1./pow(dr,12)-0.5/pow(dr,6);
	}
      }
    }
    Walker->SumComp(iV,v*4.);
    Walker->SumComp(iW,w*48./3.);
}
//Print results for current block
void CanonicEns::Averages(int iblk, ofstream *OutRes){
    double vAv,wAv,gAv,bin_m,bin_M,dVr;
    DataVett stima(n_Props);
    const int wd=12;

    if(iblk==1)	{cout << endl << endl;}
    cout << "Block number " << iblk << endl;
    cout <<"Acceptance rate "<< double(accepted)/double(tries) <<endl<<endl;
    cout << "----------------------------" << endl << endl;

    vAv=Walker->GetComp(iV)/double(nStep);
    wAv=Walker->GetComp(iW)/double(nStep);
    stima.DefComp(iV,Epot(vAv));	//Potential energy
    stima.DefComp(iW,Pressure(wAv));	//Pressure
    for(int ibin=0; ibin<nBins; ++ibin){
        gAv=Walker->GetComp(iGr+ibin)/double(nStep);
	bin_m=ibin*BinSize;		//bin's left extreme
	bin_M=bin_m+BinSize;		//bin's right extreme
	dVr=4./3.*M_PI*(pow(bin_M,3)-pow(bin_m,3));
	stima.DefComp(iGr+ibin,gAv/(Rho*nPart*dVr));//g(r)
    }

    DoAverages(stima,iblk);
    m_Props=iGr;		//write progressive res and err only for P,V
    WriteAverages(stima,iblk,OutRes);
    m_Props=n_Props;		//(modified m_Props only for WriteAverages)

    if(iblk==1)	{OutRes[iGr] << setw(wd) << "block #" << setw(2*wd) << "bins" << endl << endl;}
    OutRes[iGr]<< setw(wd) << iblk << setw(2*wd);
    for(int i=0; i<nBins; i++)	{OutRes[iGr]<< stima.GetComp(iGr+i) << setw(wd);}
    OutRes[iGr] << endl;
}


//Print simulation's information and parameters
void CanonicEns::PrintInfo() const {
    cout << endl << "Classic Lennard-Jones fluid        " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;

    cout << "Temperature = " << Temp << endl;
    cout << "Number of particles = " << nPart << endl;
    cout << "Density of particles = " << Rho << endl;
    cout << "Volume of the simulation box = " << Vol << endl;
    cout << "Edge of the simulation box = " << Box << endl;
    cout << "Cutoff of the interatomic potential = " << Rcut << endl << endl;
    cout << "Tail correction for the potential energy = " << Vtail << endl;
    cout << "Tail correction for the virial           = " << Ptail << endl;

    cout << "The program perform Metropolis moves with uniform translations" << endl;
    cout << "Moves parameter = " << Delta << endl;
    cout << "Number of blocks = " << nBlk << endl;
    cout << "Number of steps in one block = " << nStep << endl << endl;
}
