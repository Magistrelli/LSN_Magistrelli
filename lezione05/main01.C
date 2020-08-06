#include "Hydrogen.h"

int main (int argc, char *argv[]){
if (argc!=2) {
    cerr << "Usage: " << argv[0] << " <boolean far variable>" << endl;
    return -1;
}

//Variables
int M=int(1e7),N=int(1e2);	//number of throws and blocks
int L=int(M/N);			//number of throws in each block
int nburn=int(2e4);		//(unnecessary) burn-in
int nlm[3]={0,0,0};		//Hydrogen states' parameters
DataVett start(3);		//starting point
//params[0][l] width of costant T, params[1][l] sigma of Gaussian T
double params[2][2]={{1.25,3.},{0.75,2.}};//params[ig][0] for GS, [ig][1] 210
Hydrogen *Metro;		//Metropolis algorithm (dynAlloch for clean)
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");//initialization random_gen from files
ofstream OutRes[1],OutConf,OutBurn,Null[1];

bool far=atoi(argv[1]);
if(far==0){			//from the origin
    for(int j=0; j<3; j++)	{start.DefComp(j,0.);}
}
if(far==1){			//far form origin
    for(int j=0; j<3; j++)	{start.DefComp(j,1e4);}
}

//trials to find a good delta
/*
Hydrogen *Test;
double param;
int IG=0,ll=1;			//IG=0 uniform T, IG=1 Gaussian T
nlm[0]=ll+1,nlm[1]=ll,nlm[2]=0;	//ll=0 GS, ll=1 210
Test = new Hydrogen(nlm,start,rnd,IG);
Test->SetBlkSt(N,int(L/10));
cout << endl << "Test for mean acceptance <A> with 50% criterion"<<endl<<endl;
for(int i=0; i<5; i++){
    cout << "i = " << i;
    if(i!=0)	{Test->Restart(start);}
    param=(i+3.)/2.;	//(i+3.) uniform, (i+2.) Gaussian,/4 pGS, /2 p210
    Test->SetParam(param);
    Test->DoSampling(Null,Null[0],0);
    cout << endl << "param = " << param << "\t<A> = " << Test->Acceptance() << endl << endl;
}
delete Test;*/

///*
//point 1.
for(int ig=0; ig<2; ig++){		//ig=0 uniform T, ig=1 Gaussian T

  for(int l=0; l<2; l++){		//l=0 GS, l=1 210
    nlm[0]=l+1,nlm[1]=l,nlm[2]=0;
    Metro = new Hydrogen(nlm,start,rnd,ig);
    Metro->SetBlkSt(N,L);		//set #blocks, #steps per block
    Metro->SetBurn(nburn);
    Metro->SetParam(params[ig][l]);

    cout << endl << "Simulation with ";
    if(ig==0){
        cout << "Uniform T-matrix of ";
	if(far==0) OutRes[0].open("results/res_"+to_string(l)+".out");
	if(far==1) OutRes[0].open("results/res_far_"+to_string(l)+".out");
	if(l==0){if(far==0) OutConf.open("results/GS.config");
		 if(far==1) OutConf.open("results/GS_far.config");}
	else	{if(far==0) OutConf.open("results/210.config");
		 if(far==1) OutConf.open("results/210_far.config");}
    } else {
	cout << "Gaussian T-matrix of ";
	if(far==0) OutRes[0].open("results/res_gauss_"+to_string(l)+".out");
	if(far==1) OutRes[0].open("results/res_gauss_far_"+to_string(l)+".out");
    }
    if(l==0)	{cout << "Ground State";}
    if(l==1)	{cout << "(n=2, l=1, m=0) state";}
    cout << endl << endl;    

    if(far==0) OutBurn.open("results/burn_in_"+to_string(ig)+"_"+to_string(l)+".out");
    if(far==1) OutBurn.open("results/burn_in_far_"+to_string(ig)+"_"+to_string(l)+".out");
    Metro->BurnIn(OutBurn,1);
    OutBurn.close();
    Metro->DoSampling(OutRes,OutConf,(ig+1)%2); //write conf only for ig==0
    cout << endl << "Mean acceptance <A> = " << Metro->Acceptance() << endl << endl;
    cout << "----------------------------" << endl << endl;

    OutConf.close();
    OutRes[0].close();
    delete Metro;
  }

}
//*/

return 0;
}
