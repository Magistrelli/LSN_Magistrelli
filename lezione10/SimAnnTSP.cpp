#include "SimAnnTSP.h"

//ATTENTION!
//Here nPart from Metropolis and Ncit from Path are the same thing!!
//Here Gene works as Xnew for Metropolis

//constructor
MetroPath::MetroPath(Random* rnd): Metropolis(rnd) {
    ifstream ReadInput;

    ReadInput.open("input.dat");
    ReadInput >> StepMax;
    nStep=StepMax;
    Istep=0,EndStep=0;
    nBlk=1;		//not blocking method, min search
    ReadInput >> Ncit;
    nPart=Ncit;
    TotD=DimSp*Ncit;	//(DimSp=1, values are ordered cities)
    ReadInput >> Square;
    if (Square)	{Loc=new string(".square.");}
    else	{Loc=new string(".circ.");}
    *Loc+=to_string(Ncit)+".SA.";
    ReadInput >> Box;
    ReadInput.close();

    Pos=new DataVett(Ncit);
    *Pos=GenRndPos(Ncit,Square,Box,Rnd);
    PathIndex=0;	//no population
    Flag=0;
    _rnd=Rnd;		//only one random number generator;
    m_Props=n_Props;
    Gene=new int[Ncit];
    RndGen();		//starting path created randomly (init also Loss)
    MinLenght=Loss;
    InitVett();
    for(int i=0; i<Ncit; ++i) {X->DefComp(i,Gene[i]);}

    cout<<"Types of mutation's operators (city 1 excluded):"<<endl;
    cout<<"  Switch two random cities"<<endl;
    cout<<"  Shift a random number of contiguous cities a random number of positions later"<<endl;
    cout<<"  Switch rnd m contiguous cities with other m at a rnd distance (in the path)"<<endl;
    cout<<"  Invert the order of a random number of cities from a random position"<<endl;

    WriteStart(Square,Box,*Loc);
}


//Here we use every single mutation to try to change configurations
void MetroPath::Move() {
    MutSwitch();		//save new cities' sequence in Gene
    SubMove();			//check if move is accepted
    MutShift();
    SubMove();
    MutPermut();
    SubMove();
    MutInvert();
    SubMove();
    
    Istep++;
    if(Ending) {EndStep++;}
    WriteLength();
}
//check and save new config
void MetroPath::SubMove() {
    LossCost();				//measure new loss
    HNew=Loss;
    if(StdIfMove(DataVett(0),-1)){	//if accepted save new config
	for(int i=0; i<Ncit; ++i) {X->SetComp(i,Gene[i]);}
	HNow=HNew;			//save new stored energy
	if(HNew<MinLenght){
	    WriteBest(*Loc);		//overwrite best file
	    MinLenght=HNew;
	    if(Ending) {nStep=EndStep+StepMax;}//if in end part (DoSampling), proceed for other StepMax steps
	}
    } else {				//else restore old config
	for(int i=0; i<Ncit; ++i) {Gene[i]=int(X->GetComp(i));}
	Loss=HNow;
    }
}


void MetroPath::PrintAcc() const
    {cout<<"\tAcceptance: "<<double(accepted)/double(tries)<<"\tbest path loss:\t"<<setprecision(6)<<MinLenght<<endl;}


//write path's progressive lenght
void MetroPath::WriteLength() const {
    ofstream OutLength;
    const int iprint=int(1e4);
    const int wd=20;

    if(Istep==1){OutLength.open("results/path"+*Loc+"lenght");
		OutLength.close();}	//clean old output files
    OutLength.open("results/path"+*Loc+"lenght",ios::app);
    OutLength<<setw(wd)<<Istep<<setw(wd)<<Loss<<endl;
    OutLength.close();    
    if(Ending)
      if(Istep%iprint==0) {cout<<"Number of steps: "<< Istep <<"\tbest path loss:\t"<<MinLenght<<endl;}
}
