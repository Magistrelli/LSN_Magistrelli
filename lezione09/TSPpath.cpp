#include "TSPpath.h"

//constructor
Path::Path(): Ncit(0),Loss(0),PathIndex(0),Flag(0) {Pos=NULL,_rnd=NULL,Gene=NULL;}

Path::Path(DataVett* pos, Random* rnd): Ncit(0),Loss(0),PathIndex(0) {
    Pos=pos;			//equals the pointers, useless copy *pos
    _rnd=rnd;
    Gene=NULL;
    Flag=0;
}


//create starting path reshuflling randomly an ordered sequence 0...Ncit-1
//computational cost prop to N (it swap every element with a random one)
void Path::RndGen() {
    int sort;
    for(int i=0; i<Ncit; ++i) {Gene[i]=i;}//crescent order (0,1,..,N-1)
    for(int i=1; i<Ncit; ++i){
	sort=int(_rnd->Rannyu(1.,Ncit));
	SwitchGene(i,sort);
    }
    Check();			//check if not errors
    LossCost();			//save chromosome's loss function
}


//Switch of two random cities
void Path::MutSwitch() {
    int city1=int(_rnd->Rannyu(1.,Ncit));	//random selected city
    int city2=int(_rnd->Rannyu(1.,Ncit));	//the 1-st is excluded
    SwitchGene(city1,city2);
    Check();
}

//shift rnd #m of contiguous cities at rnd starting point ini at rnd dist n
//0 must be fixed, so min and max value of m are [2,Ncit-2] (Ncit-1 would be all fixed)
//rnd() is in [0,1), smaller p ->  bigger m -> operator less destructive
//the m cities must be contiguous, so max value of ini is Ncit-m (no PBC)
//fixed 0 destroy contiguity of the m cities if cicle, so not use PBC
//if n<0 shitf at left, if n>0 shift at right, if n=0 sort again
void Path::MutShift(){
    int m,ini,n;
    double const p=0.5;
    bool right=1;		//if true shift at right, else at left

    m=2+int((Ncit-3.)*pow(_rnd->Rannyu(),p));
    ini=int(_rnd->Rannyu(1.,Ncit-m+1));
    do {n=int(_rnd->Rannyu(-ini,Ncit-(ini+m)+1));}//e.g. int(-2.5)=-2
    while(n==0);
    if(n<0) {right=0;}
//  cout<<endl<<"ini "<<ini<<" m "<<m<<" n "<<n<<endl; //usefull for check

    if(right) {
	for(int j=0; j<n; ++j){
	    for(int i=ini+m; i>ini; --i) {SwitchGene(i-1,i);}
	    ini++;
	}
    } else {
	for(int j=0; j>n; --j){
	    for(int i=ini; i<ini+m; ++i) {SwitchGene(i-1,i);}
	    ini--;
	}
    }
    
    Check();
}

//switch rnd #m of contiguous cities at rnd pos ini with other m at rnd pos n
//rnd() is in [0,1), smaller p ->  bigger m -> operator less destructive
void Path::MutPermut(){
    int m,ini,n;
    double const p=0.5;

    m=2+int((Ncit/2.-2)*pow(_rnd->Rannyu(),p));//need to switch, m < Ncit/2
    do {ini=int(_rnd->Rannyu(1.,Ncit-m+1));}//same thoughts as in MutShift
    while((ini-1<m)&&(Ncit-(ini+m)<m));	//not enough contig cities to switch
    do {n=int(_rnd->Rannyu(1,Ncit-m+1));}//max pos has n+m=Ncit
    while ((n>ini-m)&&(n<ini+m));	//other different m cities
//  cout<<endl<<"ini "<<ini<<" m "<<m<<" n "<<n<<endl; //usefull for check

    for(int i=0; i<m; ++i) {SwitchGene(ini+i,n+i);}
    Check();
}

//invert order of cities from rnd position ini to rnd pos end
void Path::MutInvert(){
    int ini=int(_rnd->Rannyu(1,Ncit-1));		//max is Ncit-2
    int end=int(_rnd->Rannyu(ini+1,Ncit));	//at least 2 cities
//  cout<<endl<<"ini "<<ini<<" end "<<end;	//usefull for check
    for(int i=0; i<=int((end-ini)/2.); ++i) {SwitchGene(ini+i,end-i);}
    Check();
}

//switch two genes' values
void Path::SwitchGene(int c1, int c2) {
    int appo=Gene[c1];
    Gene[c1]=Gene[c2];
    Gene[c2]=appo;
}


//cost/loss function evaluation
void Path::LossCost() {
    double loss=0.;
    for(int i=0; i<Ncit-1; ++i)	{loss+=Distance(Gene[i],Gene[i+1]);}
    loss+=Distance(Gene[Ncit-1],Gene[0]); //dist from last city and first one
    Loss=loss;
}
//distance between two cities
double Path::Distance(int i1, int i2) const{
    double dist=0.;
    for(int j=0; j<Nspace; ++j) {dist+=pow(Pos->GetComp(Nspace*i1+j)-Pos->GetComp(Nspace*i2+j),2);}
    return sqrt(dist);
}


//check if every city appears one and only one time
void Path::Check() {
    int city;			//encountered city
    bool check[Ncit];
    for(int i=0; i<Ncit; ++i) {check[i]=0;}

    if(Gene[0]!=0) {
	cout<<endl<<"ERROR! The "<<PathIndex+1<<"-th chromosome hasn't 1 as first city!"<<endl;
	exit(-2);
    }
    for(int i=0; i<Ncit; ++i) {	//over all genes
      city=Gene[i];		//city 1 is 0, 2 is 1, etc
      if(check[city]) {		//city i encountered yet
	cout<<endl<<"ERROR! city "<< city+1 <<" repeated!" << endl;
	cout<<"  The "<<PathIndex+1<<"-th chromosome will be deleted and substituted by a randomly generated new one"<<endl;

	RndGen();
	Flag=1;			//error encountered
	i=Ncit;			//useless continue with check
      } else	{check[city]++;}//city encountered for the first time
    }
}


//write starting path in default file
void Path::WriteStart(bool square, double box, string loc) const {
    ofstream OutPos,OutPath;
    cout<<endl<<"Cities' position are randomly generated ";
    if(square)	{cout<<"in a square of side = " << 2.*box << endl;}
    else	{cout<<"on a circumference of radius = "<<box<<endl;}
    cout<<"Saved cities' positions in file results/cities"+loc+"pos"<<endl;

    OutPos.open("results/cities"+loc+"pos");
    for(int i=0; i<Ncit; ++i){
	for(int j=0; j<Nspace; ++j) {OutPos << Pos->GetComp(Nspace*i+j) << "   ";}
	OutPos << endl;
    }
    OutPos.close();

    cout<<"Saved starting path in file results/cities"+loc+"start"<<endl;
    OutPath.open("results/cities"+loc+"start");
    WritePath(OutPath);
    OutPos.close();
}
//write best path in default file
void Path::WriteBest(string loc) const {
    ofstream OutPath;
    OutPath.open("results/path"+loc+"best");
    WritePath(OutPath);
    OutPath.close();
}
//write path of 1-st chromosome
void Path::WritePath(ofstream& OutPath) const {
    OutPath<<Loss<<endl<<endl;
    for(int i=0; i<Ncit; ++i){
	for(int j=0; j<Nspace; ++j) {OutPath << Pos->GetComp(Nspace*Gene[i]+j) << "   ";}
	OutPath << endl;
    }
    OutPath << Pos->GetComp(Nspace*Gene[0]) <<"   "<< Pos->GetComp(Nspace*Gene[0]+1);	//end in first city
}





//EXTERNAL FUNCTION

//generate cities' random positions
DataVett GenRndPos(int ncit, bool square, double box, Random* rnd) {
    int totD=ncit*Nspace;
    DataVett pos(totD);		//(1x,1y,2x,2y,...,Ncit_x,Ncit_y)
    
    cout<<"Compute cities' positions randomly" << endl;
    if(square) {		//cities in a square
	for(int i=0; i<totD; ++i) {pos.DefComp(i,rnd->Rannyu(-box,box));}
    } else {			//cities on a circumference
	for(int i=0; i<ncit; ++i) {
	    double theta=rnd->Rannyu(0.,2.*M_PI);//random angle in rad
	    pos.DefComp(Nspace*i,box*cos(theta));	//ix
	    pos.DefComp(Nspace*i+1,box*sin(theta));	//iy
	}
    }
    return pos;
}
