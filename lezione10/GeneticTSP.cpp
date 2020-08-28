#include "GeneticTSP.h"

//CHROMO

//constructor
Chromo::Chromo(): Path(),MutProb(0) {_rnd=NULL;}
Chromo::Chromo(DataVett* pos, Random* rnd, int ichr, int rank): Path(pos,rnd) {
    double appo;	//variable to discard useless input
    ifstream ReadInput;
    
    ReadInput.open("input.dat");
    ReadInput >> appo;	//GenMax already in Population
    ReadInput >> appo;	//Nmigr already in Population
    ReadInput >> Ncit;
    ReadInput >> appo;	//Nchr already in Population
    ReadInput >> appo;	//Square already in Population
    ReadInput >> appo;	//Box already in Population
    ReadInput >> appo;	//CrossProb already in Population
    MutProb=new double[NMut];
    for(int i=0; i<NMut; ++i) {ReadInput >> MutProb[i];}
    ReadInput.close();

    PathIndex=ichr;	//chromosome's index in population
    Gene=new int[Ncit];	//Ncit genes in a chromosomes, int values
    RndGen();		//starting path created randomly (init also Loss)

    if((PathIndex==0)&&(rank==0)) {//rank for parallel computing (given from population)
      cout<<"Compute starting chromosomes randomly" << endl;
      cout<<"Types of mutation's operators (city 1 excluded):"<<endl;
      cout<<"  Switch two random cities, probability = "<< MutProb[0] << endl;
      cout<<"  Shift a random number of contiguous cities a random number of positions later, prob = " << MutProb[1] << endl;
      cout<<"  Switch rnd m contiguous cities with other m at a rnd distance (in the chromo), prob = " << MutProb[2] << endl;
      cout<<"  Invert the order of a random number of cities from a random position, prob = " << MutProb[3] << endl;
    }
}

//define operator = for Chromo
Chromo& Chromo::operator=(const Chromo& chr){
    Ncit=chr.Ncit;
    Pos=chr.Pos;
    _rnd=chr._rnd;
    PathIndex=chr.PathIndex;
    Loss=chr.Loss;
    Flag=chr.Flag;

    if(MutProb) {delete[] MutProb;}
    MutProb=new double[NMut];
    for(int i=0; i<NMut; ++i) {MutProb[i]=chr.MutProb[i];}
    if(Gene) {delete[] Gene;}
    Gene=new int[Ncit];
    for(int i=0; i<Ncit; ++i) {Gene[i]=chr.Gene[i];}
    
    return *this;
}

//try all possible mutation
void Chromo::Mutation(){
    double sort;
    sort=_rnd->Rannyu();
    if(sort<MutProb[0]) {MutSwitch();}
    sort=_rnd->Rannyu();
    if(sort<MutProb[1]) {MutShift();}
    sort=_rnd->Rannyu();
    if(sort<MutProb[2]) {MutPermut();}
    sort=_rnd->Rannyu();
    if(sort<MutProb[3]) {MutInvert();}
}






//POPULATION

//constructor
Population::Population() {
    ifstream ReadInput;
    ReadInput.open("input.dat");
    ReadInput >> GenMax;
    ReadInput >> Nmigr;
    ReadInput >> Ncit;
    TotD=Ncit*Nspace;
    ReadInput >> Nchr;
    ReadInput >> Square;
    if (Square)	{Loc=new string(".square.");}
    else	{Loc=new string(".circ.");}
    *Loc+=to_string(Ncit);
    ReadInput >> Box;
    ReadInput >> CrossProb;
    ReadInput.close();
    ParSize=1,ParRank=0;	//no parallelized by default
    MinLoss=1e3;		//very big value
}
Population::Population(Random* rnd): Population() {
    Rnd=rnd;
    *Loc+=".";
    CommConstr();
}
Population::Population(Random* rnd, int size, int rank): Population() {
    Rnd=rnd;
    ParSize=size,ParRank=rank;
    *Loc+="_"+to_string(rank)+".";
    CommConstr();
}
//constructors' common operations
void Population::CommConstr() {
    double *poscit;		//usefull for cores communication
    
    Pos=new DataVett(TotD);
    if(ParRank==0){
	cout<<endl<< "Travelling Salesman Problem" << endl;
	cout<< "Genetic Algorithm simulation" << endl << endl;
	cout<< "Max generations waiting for improving    = "<< GenMax <<endl;
	if(ParSize>1) cout<<"Generations between migrations           = "<<Nmigr<<endl;
	cout<< "Number of Cities                         = "<< Ncit <<endl;
	cout<< "Number of Chromosomes in each Generation = "<< Nchr <<endl;
	*Pos=GenRndPos(Ncit,Square,Box,Rnd);
	cout<<endl<<"crossover remove final part of 2 chromosome and refill with genes' order of the other one with probability "<<CrossProb<<endl;
    }
    if(ParSize>1){
	poscit=new double[TotD];
	if(ParRank==0) for(int i=0; i<TotD; ++i) {poscit[i]=Pos->GetComp(i);}
	MPI_Bcast(&poscit[0],TotD,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD);
	if(ParRank!=0) for(int i=0; i<TotD; ++i) {Pos->DefComp(i,poscit[i]);}
	delete[] poscit;
    }
    
    Chr=new Chromo*[Nchr];	//initialize chromosomes
    for(int i=0; i<Nchr; ++i) {Chr[i]= new Chromo(Pos,Rnd,i,ParRank);}

    Chr[0]->WriteStart(Square,Box,*Loc);
}
//destructor
Population::~Population() {
    delete[] Chr;
    delete Loc;
    delete Pos;
}


//Simulate genetic evolution
void Population::Evolution() {
    int sel[2],igen=1,nnmig=0;	//parents chromosomes number, generation index, migration counter
    int count=0,incount=0;	//#gen of non improving MinLoss
    Chromo newChr[Nchr];	//brand new generation
    double *mins,writtenMin=1e3;//min recv (parallel code), min in this core's file best

    if(Chr[0]->GetLoss()<MinLoss) {MinLoss=Chr[0]->GetLoss();}
    if(ParRank==0) cout<<endl<<"init 1st path loss:\t\t\t\t"<<MinLoss<<endl;


    while(count<GenMax){	//count is modified only one gen before a migration

      for(int i=0; i<Nchr; i+=2){	//sons generated in pairs
	for(int j=0; j<2; ++j) {sel[j]=Selection();}//parents' selection
	Crossover(sel[0],sel[1]);	//sons S0, S1 generated by crossover (maybe)
	for(int j=0; j<2; ++j) {S[j].Mutation();}//try all mutations on sons
	newChr[i]=S[0];			//construct new generation
	if(i==Nchr-1){}			//if we are here Nchr is odd, discard last son
	else {newChr[i+1]=S[1];}
      }
      for(int i=0; i<Nchr;++i){
	*Chr[i]=newChr[i];	//replace old generation
	Chr[i]->LossCost();	//and evaluate new loss function's values
      }
      Ordering();		//1st chromosome will be the best one
      if((ParSize>1)&&(igen%Nmigr==0)){
	  Migration(nnmig);
	  if(Chr[0]->GetLoss() > Chr[1]->GetLoss()) {Ordering();}
	  nnmig++;
      }

      if(Chr[0]->GetLoss()<MinLoss){	//new best path finded
	Chr[0]->WriteBest(*Loc);	//overwrite best file
	MinLoss=Chr[0]->GetLoss();
	if(ParSize>1) writtenMin=MinLoss;
	count=0,incount=0;
      } else {incount++;}
      if((ParSize>1)&&(Chr[0]->GetLoss()<writtenMin)){//maybe newmin>=MinLoss but newmin<writtenMin, so change also this core's best file
	Chr[0]->WriteBest(*Loc);	//here Chr[0] is that from migration
	writtenMin=Chr[0]->GetLoss();
      }

      if((ParSize>1)&&((igen+1)%Nmigr==0)){//one gen before a migration, see if someone else has a better MinLoss
        mins=new double[ParSize];
	MPI_Allgather(&MinLoss,1,MPI_DOUBLE_PRECISION,&mins[0],1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD);
	for(int ic=0; ic<ParSize; ++ic){
	  if(mins[ic]<MinLoss){	//not for the core that found the min
	    MinLoss=mins[ic];
	    count=0,incount=0;
	  }
	}
	delete[] mins;
      }
      if((igen+1)%Nmigr==0){	//one gen before Migration
	count+=incount;		//update count only here
	incount=0;
      }

      WriteLength(igen);	//update instant lenght of best path and 1st half best path, cout progression
      igen++;
    }

    
    cout<<endl<<"Saved 1st and 1st half paths lengths evolution in file results/path"+*Loc+"lenght"<<endl;
    cout<<"Saved best path in file results/path"+*Loc+"best"<<endl;
    WriteEnd();			//write final configuration (all chromosomes)
}

//Simulate genetic evolution
void Population::Migration(int nnmig) {
    MPI_Status stat1,stat2,stat3,stat4;
    MPI_Request req1,req2,req3,req4;
    int *mesg1=new int[Ncit],*mesg2=new int[Ncit];	//msg is 1st chromosome of the continents
    double lMsg1,lMsg2;					//and its lenght
    int tag[4];
    for(int i=0; i<4; ++i) {tag[i]=ParSize+4*nnmig+i;}	//different from ic of Evol and different every time
    int core[ParSize],sort,appo;

    if(ParRank==0){
      for(int i=0; i<ParSize; ++i) {core[i]=i;}	//crescent order
      for(int i=0; i<ParSize; ++i){		//random permutation for random migration
	sort=int(Rnd->Rannyu(0.,ParSize));
	appo=core[i];
	core[i]=core[sort];
	core[sort]=appo;
      }
    }
    MPI_Bcast(core,ParSize,MPI_INTEGER,0,MPI_COMM_WORLD);
    
    for(int ir=1; ir<ParSize; ir+=2){	//run over couple of cores, not last rank if ParSize is odd
      if(ParRank==core[ir]){
	for(int ic=0; ic<Ncit; ++ic) {mesg1[ic]=Chr[0]->GetGene(ic);}
	MPI_Isend(&mesg1[0],Ncit,MPI_INTEGER,core[ir-1],tag[0],MPI_COMM_WORLD,&req1);
	MPI_Irecv(&mesg2[0],Ncit,MPI_INTEGER,core[ir-1],tag[1],MPI_COMM_WORLD,&req2);
	lMsg1=Chr[0]->GetLoss();
	MPI_Isend(&lMsg1,1,MPI_DOUBLE_PRECISION,core[ir-1],tag[2],MPI_COMM_WORLD,&req3);
	MPI_Irecv(&lMsg2,1,MPI_DOUBLE_PRECISION,core[ir-1],tag[3],MPI_COMM_WORLD,&req4);

	MPI_Wait(&req2,&stat2);		//wait new 1st chromo (all communications here are non-blocking)
	for(int ic=0; ic<Ncit; ++ic) {Chr[0]->SetGene(ic,mesg2[ic]);}
	MPI_Wait(&req4,&stat4);		//and its lenght
	Chr[0]->SetLoss(lMsg2);		//overwrite it
      } 
      else if(ParRank==core[ir-1]){
	MPI_Irecv(&mesg1[0],Ncit,MPI_INTEGER,core[ir],tag[0],MPI_COMM_WORLD,&req1);
	for(int ic=0; ic<Ncit; ++ic) {mesg2[ic]=Chr[0]->GetGene(ic);}
	MPI_Isend(&mesg2[0],Ncit,MPI_INTEGER,core[ir],tag[1],MPI_COMM_WORLD,&req2);
	MPI_Irecv(&lMsg1,1,MPI_DOUBLE_PRECISION,core[ir],tag[2],MPI_COMM_WORLD,&req3);
	lMsg2=Chr[0]->GetLoss();
	MPI_Isend(&lMsg2,1,MPI_DOUBLE_PRECISION,core[ir],tag[3],MPI_COMM_WORLD,&req4);

	MPI_Wait(&req1,&stat1);		//wait new 1st chromo (here all ir-1's Isend are already done)
	for(int ic=0; ic<Ncit; ++ic) {Chr[0]->SetGene(ic,mesg1[ic]);}
	MPI_Wait(&req3,&stat3);		//and its lenght
	Chr[0]->SetLoss(lMsg1);
      }
    }
}


//crossover operator
//for each of the two chromosomes delete all the cities from rnd position pos and rewrite them in the same order of the other chromosome
//pos is in [2,Ncit-2]: pos=1 would only exchange PathIndex, pos=N_cit-1 would do nothing
void Population::Crossover(int sel1, int sel2) {
    bool equal=1;
    int pos=int(Rnd->Rannyu(2.,Ncit-1.));//from 2 to N_cit-2
    int m=Ncit-pos,search,find;		//#cities deleted, running index, #found cities in the other chromosome
    int appo[2*m];
    Chromo* p[2]={Chr[sel1],Chr[sel2]};	//parents
    S[0]=*p[0],S[1]=*p[1];		//init sons as copy of parents


    if(Rnd->Rannyu()<CrossProb){//if not will return parents again
      if(sel1!=sel2){		//check if parents are identical
	for(int i=1; i<Ncit; ++i) {
	    if(p[0]->GetGene(i)!=p[1]->GetGene(i)){
		equal=0;	//parents are different
		i=Ncit;		//stop for cicle
	    }
	}
      }
      if(equal==0){		//if p1==p2 do nothing

	for(int i=0; i<2; ++i){ //both the chromosome
	  search=1,find=0;
	  while(find<m){	//we have to find all the m cities in the other chromosome

	    if(search>=Ncit){
		cerr<<endl<<"Error in crossover! removed elements from "<<i+1<<" chromo not found in the other one!!"<<endl;
		exit(-3);
	    }
	    
	    for(int j=0; j<m; ++j){//see if other chromo's search-th gene is one of the m cities of this chromo (for both, i=0 and i=1)
		if(p[(i+1)%2]->GetGene(search)==p[i]->GetGene(pos+j)){//for i=0 search in the 2nd chromosome
		    appo[i*m+find]=p[i]->GetGene(pos+j);//save the new order
		    find++;	//found one of the m cities
		    j=m;	//exit for cicle
		}
	    }
	    search++;		//next other chromo's element
	  }
	  for(int j=0; j<m; ++j) {S[i].SetGene(pos+j,appo[i*m+j]);}
	  S[i].Check();
	}

      }
    }
}


//Select individual with weight 1/e^(loss) (ordering not needed)
//the better is a path respect to others the more it will reproduce itself
int Population::Selection() const {
    double fit[Nchr],sort,sum=0.;//fitness function, sorted number
    for(int i=0; i<Nchr; ++i){
	sum+=1./exp(Chr[i]->GetLoss());//can be changed (power law, ...)
	fit[i]=sum;		//accumulated probability (not norm)
    }
    sort=Rnd->Rannyu(0.,sum);
    for(int i=0; i<Nchr; ++i) {
	if(sort<fit[i]) {return i;}
    }
    return -1;			//if this happens ther's an error
}


//fitness based population ordering, path with less cost come first
//I use BubbleSort, not the fastest but simple and elements to order are few
void Population::Ordering() {
    bool exc=1;
    Chromo appo;
    
    for(int i=0; i<Nchr-1; ++i) {
      while(exc){	//if not switch in an iteration then ordering done
      exc=0;
	for(int j=0; j<Nchr-i-1; ++j) {	//last i elements already in place
	  if(Chr[j]->GetLoss() > Chr[j+1]->GetLoss()) {
	    appo=*Chr[j];
	    *Chr[j]=*Chr[j+1];
	    *Chr[j+1]=appo;
	    exc=1;
	  }
	}
      }
    }
}


//write lenght of 1-st path and mean of 1-st half
void Population::WriteLength(int igen) const {
    double sum=0.;
    ofstream OutLoss;
    const int iprint=GenMax/10;
    const int wd=20;

    for(int i=0; i<int(Nchr/2.); ++i) {sum+=Chr[i]->GetLoss();}
    if(igen==1){OutLoss.open("results/path"+*Loc+"lenght");	//delete file from previous simulation
		OutLoss.close();}
    OutLoss.open("results/path"+*Loc+"lenght",ios::app);
    OutLoss<<setw(wd)<<igen<<setw(12)<<Chr[0]->GetLoss()<<setw(wd)<<sum/int(Nchr/2.)<<endl;
    OutLoss.close();    
    if((ParRank==0)&&(igen%iprint==0)) {cout<<"Number of generations: "<< igen <<"\tbest path loss:\t"<<MinLoss<<endl;}
}

//write final chromosomes and final seed in default files
void Population::WriteEnd() {
    ofstream OutChromo;
    cout<<"Saved final chromosomes in file results/chromo"+*Loc+"end"<<endl;
    OutChromo.open("results/chromo"+*Loc+"end");
    for(int i=0; i<Nchr; ++i){
	for (int j=0; j<Ncit; ++j) {OutChromo<<setw(12)<<Chr[i]->GetGene(j);}
	OutChromo<<endl;
    }
    OutChromo.close();
    Rnd->SaveSeed();
}
