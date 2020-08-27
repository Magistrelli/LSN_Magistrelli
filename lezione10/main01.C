#include "SimAnnTSP.h"

int main (int argc, char *argv[]){
cout<<endl;

//Variables
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");	//init random_gen from files
MetroPath path(rnd);
const int nphases=5, ncit=path.GetNcit();//#phases of cooling, #cities
int nbeta[nphases];			//#different beta in each phase
int nbstep[nphases];			//#step for each beta
double Tmin[nphases];			//min T for each phase
if (ncit==32){
    nbeta[0]=4,nbeta[1]=6,nbeta[2]=6,nbeta[3]=5,nbeta[4]=10; //#annealing_beta
    nbstep[0]=300,nbstep[1]=600,nbstep[2]=1000,nbstep[3]=1500,nbstep[4]=3000;
    Tmin[0]=2.5,Tmin[1]=1.,Tmin[2]=0.3,Tmin[3]=0.1,Tmin[4]=0.01;
    }
else if (ncit==100){
    nbeta[0]=4,nbeta[1]=6,nbeta[2]=6,nbeta[3]=10,nbeta[4]=20;
    nbstep[0]=1000,nbstep[1]=2000,nbstep[2]=3000,nbstep[3]=5000,nbstep[4]=10000;
    Tmin[0]=2.5,Tmin[1]=1.,Tmin[2]=0.3,Tmin[3]=0.1,Tmin[4]=0.001;
    }
else{
    cerr << "The program works only with Ncit=32 or Ncit=100!" << endl;
    return -1;
}
double Tmax=30.;				//max annealing T
int jj,totNbeta=0;
for(int i=0; i<nphases; ++i) {totNbeta+=nbeta[i];}
DataVett beta(totNbeta+1),bStep(totNbeta+1);	//annealing schedule

//simulation
totNbeta=0;					//betas done
for(int iphase=0; iphase<nphases; ++iphase){	//T-unif phases
    for(int i=0; i<nbeta[iphase]; ++i){
        jj=totNbeta+i;
	beta.DefComp(jj,1./(Tmax-i*(Tmax-Tmin[iphase])/nbeta[iphase]));
	bStep.DefComp(jj,nbstep[iphase]);
    }
    totNbeta+=nbeta[iphase];
    Tmax=Tmin[iphase];				//new Tmax
}
beta.DefComp(totNbeta,1./Tmin[nphases-1]);
bStep.DefComp(totNbeta,nbstep[nphases-1]);

cout<<endl;
path.SimAnnealing(beta,bStep);
cout<<endl;

delete rnd;
return 0;
}
