#include "Vettore.h"
#include "random.h"
#include <iomanip>

int main (int argc, char *argv[]){
//Variables
int M=int(1e4);		//number of means measured
int N[]={0,1,2,10,100};	//number of draws mediated
int L=N[4]*M;		//number of total extractions needed
double lambda=1.;	//exp distribution parameter
double mu=0.,Gamma=1.;	//Cauchy-Lorentz distribution parameters
double sum;
int count;		//will count numbers of draws from each distribution
ofstream output;
Random rnd;
rnd.SetRandom("Primes","seed.in");//initialization random_gen from files

//Point 1. in random.h and %.cpp, definition of rnd.Exp and rnd.Lorentz 
//Point 2.
output.open("res02.out");
for(int distr=0; distr<3; distr++){
    count=0;
    while(count<L){
	sum=0.;
	for(int j=1; j<5; j++){	//at step j will be summed in total N[j] draws
	    for(int i=N[j-1]; i<N[j]; i++){
		if(distr==0)	 {sum+=rnd.Rannyu();}
		else if(distr==1){sum+=rnd.Exp(lambda);}
		else if(distr==2){sum+=rnd.Lorentz(mu,Gamma);}
		count++;
	    }
	    output << sum/N[j] << setw(15);//S(N[j]) all in the same line
	}
	output << endl;	//endl for the next measures of S(N[0]),...,S(N[3])
    }
    output << endl << endl;
}
output.close();
return 0;
}
