#include "Ising.h"

#define nSpin nPart	//number of spin = number of degrees of freedom
#define S X		//configuration

//constructors
Ising::Ising(Random* rnd): Metropolis(rnd){//Prepare all stuff for simulation
    ifstream ReadInput,ReadConf;
    double in;

    m_Props=n_Props;
    iPrint=int(1e5);
    ReadInput.open("input.dat");
    ReadInput >> Temp;
    Beta = 1./Temp;
    ReadInput >> nSpin;
    TotD=nSpin;			//DimSp=1 (default setting)
    ReadInput >> J;
    ReadInput >> h;
    ReadInput >> UseMetro; 	//if=1 Metropolis else Gibbs
    ReadInput >> nBlk;
    ReadInput >> nStep;
    ReadInput >> _restart;

    InitVett();
    if(_restart==1){
      ReadConf.open("config.0");
      if (ReadConf.fail()) {
        cout << endl <<"File config.0 doesn't exist!"<< endl <<"\tSimulation will create initial configuration randomly (infinite T)"<< endl;
        _restart=0;
      }
      else{
	cout << "Read starting configuration from file config.0 " << endl << endl;
	for (int i=0; i<nSpin; ++i){
	    ReadConf >> in;
	    S->DefComp(i,in);
	}
      }
      ReadConf.close();
    }
    if(_restart==0){
      for (int i=0; i<nSpin; ++i){	//initial configuration, infinite T
	if(Rnd->Rannyu()>=0.5)	{S->DefComp(i,1);}
	else 			{S->DefComp(i,-1);}
      }
    }
    ReadInput.close();

    Measure();	//Evaluate energy etc. of the initial configuration
    cout << "Initial energy = " << Walker->GetComp(iU)/double(nSpin) << endl;
}


//one step of Markov Chain
void Ising::Move(){
    int o;
    bool ifMove;
    for(int i=0; i<nSpin; ++i){
	o=int(Rnd->Rannyu()*nSpin);	//Select spin randomly, 0<=o<=nSpin-1
	if(UseMetro==1){		//Metropolis sampling
	    ifMove=StdIfMove(*X,o);	//doesn't matter first argument
	    if (ifMove==1) {S->SetComp(o,-1.*S->GetComp(o));}//spin o flipped
	} else	{Gibbs(o);}		//Gibbs sampling
    }
}
//T(x|y)=T(y|x) => q=p(Xnew)/p(Xold)  (here Snew is useless)
double Ising::qRatio(const DataVett Snew, int ip) const{
    double sm=S->GetComp(ip);	//ip=iflip
    double diffE=2.*sm*(J*(S->GetComp(Pbc(ip-1))+S->GetComp(Pbc(ip+1)))+h);
    return exp(-Beta*diffE);
}

//one step of Markov Chain, Gibbs algorithm
void Ising::Gibbs(int o){
    double sumS=S->GetComp(Pbc(o-1))+S->GetComp(Pbc(o+1));
    double prob=1./(1.+exp(-2.*Beta*(J*sumS+h)));
    double sort=Rnd->Rannyu();
    if(sort<=prob)	{S->SetComp(o,1.);}
    else		{S->SetComp(o,-1.);}
    accepted++;		//Acceptance is always A=1
    tries++;
}


//Evaluate energy etc. of the actual configuration
void Ising::Measure(){
    double u=0.,m=0.;
    for (int i=0; i<nSpin; ++i){	//cycle over spins
	u += -J*S->GetComp(i)*S->GetComp(Pbc(i+1)) - h/2.*(S->GetComp(i)+S->GetComp(Pbc(i+1)));
	m+=S->GetComp(i);
    }
    Walker->SumComp(iU,u);
    Walker->SumComp(iC,u*u);
    Walker->SumComp(iM,m);
    Walker->SumComp(iX,m*m);
}
//Print results for current block
void Ising::Averages(int iblk, ofstream *OutRes){
    double eneAv,eneAv2,spinAv,spinAv2,a;
    DataVett stima(n_Props);

    eneAv=Walker->GetComp(iU)/double(nStep);
    eneAv2=Walker->GetComp(iC)/double(nStep);
    spinAv=Walker->GetComp(iM)/double(nStep);
    spinAv2=Walker->GetComp(iX)/double(nStep);

    a=eneAv/double(nSpin);
    stima.DefComp(iU,a);
    stima.DefComp(iC,Beta*Beta*(eneAv2-eneAv*eneAv)/double(nSpin));
    a=spinAv/double(nSpin);
    stima.DefComp(iM,a);
    stima.DefComp(iX,Beta*(spinAv2/double(nSpin)-a*a));

    DoAverages(stima,iblk);
}


//Algorithm for periodic boundary conditions
int Ising::Pbc(int i) const{
    if(i >= nSpin) i -= nSpin;
    else if(i < 0) i += nSpin;
    return i;
}


//Print simulation's information and parameters
void Ising::PrintInfo() const{
    cout << endl << "Classic 1D Ising model     " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;
    
    cout << "Temperature = " << Temp << endl;
    cout << "Number of spins = " << nSpin << endl;
    cout << "Exchange interaction = " << J << endl;
    cout << "External field = " << h << endl << endl;
    
    if(UseMetro==1){cout << "The program perform Metropolis moves" << endl;}
    else	   {cout << "The program perform Gibbs moves" << endl;}
    cout << "Number of blocks = " << nBlk << endl;
    cout << "Number of steps in one block = " << nStep << endl << endl;
}

#undef nSpin
#undef S
