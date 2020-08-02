#include "Vettore.h"
#include "random.h"


int main (int argc, char *argv[]){

//Variables
int M=int(1e6),N=int(1e2);	//number of throws and blocks (below intervals)
int L=int(M/N);			//number of throws in each block
DataVett mean(N),res(N),err(N);	//single mean measures, results and errors
DataVett sigma(N),resS(N),errS(N);//DataVett for the standard deviation
double sum,sumS,sort;		//additional variables
int count;
DataVett Chi(N);		//Chi^2
double exp=L/N;			//expected value in an interval
double min_int=0.,max_int=1./N;	//min and max of the interval
ofstream output;

Random rnd;
rnd.SetRandom("Primes","seed.in");//initialization random_gen from files

Vett Sort(M);
for(int i=0; i<M; i++)	{Sort.DefComp(i,rnd.Rannyu());}


//Point 1. and 2.
for(int i=0; i<N; i++){
    sum=0.;
    sumS=0.;
    for(int j=0; j<L; j++){
    	sum+=Sort.GetComp(j+i*L);
    	sumS+=pow(Sort.GetComp(j+i*L)-0.5,2);
    }
    mean.DefComp(i,sum/L);
    sigma.DefComp(i,sumS/L);
}
mean.StatErrProg(res,err);	//computation of res ed err, blocking method
sigma.StatErrProg(resS,errS);

//Point 3.
for(int i=0; i<N; i++){	//cycle on intervals
    count=0;
    for(int j=0; j<L; j++){
	sort=Sort.GetComp(j+i*L);	//for interval i use numbers in block i
	if(min_int<sort && sort<max_int)	{count++;}
    }
    if(i==0)	{Chi.DefComp(i,pow(count-exp,2)/exp);}
    else	{Chi.DefComp(i,Chi.GetComp(i-1)+pow(count-exp,2)/exp);}
    min_int=max_int;	//go to the next interval
    max_int+=1./N;
}

//writing results
output.open("res01.out");
output << M << "," << N << endl << endl;
for (int i=0; i<N; i++)	{output << res.GetComp(i) << "," << err.GetComp(i) << "," << resS.GetComp(i) << "," << errS.GetComp(i) << "," << Chi.GetComp(i) << endl;}
output.close();

return 0;
}
