#include "Vettore.h"
#include "random.h"

int main (int argc, char *argv[]){

//Variables
int M=int(1e5),N=int(1e2);	//number of throws and blocks
int L=int(M/N);			//number of throws in each block
DataVett Int(N),res(N),err(N);	//single integral measures, results and errors
DataVett Imp(N),resIm(N),erIm(N);//DataVett for importance sampling
double x,y,sum,suIm;		//additional variables
double dmax=M_PI/2.*(1./2+2./M_PI);//max value of distribution d(x)
ofstream output;

Random rnd;
rnd.SetRandom("Primes","seed.in");//initialization random_gen from files

//Point 1. and 2.
for(int i=0; i<N; i++){
    sum=0.;			//progressive sum of g(x)
    suIm=0.;			//progressive sum of g(x)/d(x)
    for(int j=0; j<L; j++){
    	x=rnd.Rannyu();		//sorted number from uniform distribution
    	sum+=M_PI/2.*cos(M_PI/2.*x);
	do {x=rnd.Rannyu();	//sorted number from d(x) with accept-reject
	    y=rnd.Rannyu(0.,dmax);}
	while (y>=M_PI/2.*(1./2+2./M_PI-x));
	suIm+=cos(M_PI/2.*x)/(1./2+2./M_PI-x);
    }
    Int.DefComp(i,sum/L);
    Imp.DefComp(i,suIm/L);
}
Int.StatErrProg(res,err);	//computation of res ed err, blocking method
Imp.StatErrProg(resIm,erIm);

//writing results
output.open("res01.out");
output << N << endl << endl;
for (int i=0; i<N; i++)	{output << res.GetComp(i) << "," << err.GetComp(i) << "," << resIm.GetComp(i) << "," << erIm.GetComp(i) << endl;}
output.close();

return 0;
}
