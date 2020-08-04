#include "Vettore.h"
#include "random.h"

int main (int argc, char *argv[]){

//Variables
int M=int(1e6),N=int(2.5e2);	//number of exstimates and blocks
int L=int(M/N);			//number of exstimates in each block
double S0=100.,K=100.;		//asset price at t=0 and strike price
double T=1.;			//delivery time
double r=0.1,sigma=0.25;	//risk-free interest rate, volatility
double Nt=int(1e2);		//number of time sub-intervals
DataVett C(N),P(N);		//single price measures, call and put
DataVett resC(N),erC(N);	//single results and errors, call
DataVett resP(N),erP(N);	//single results and errors, put
DataVett resC_d(N),erC_d(N);	//discretized case
DataVett resP_d(N),erP_d(N);	//discretized case
double WT,ST;			//Wiener process and asset price at t=T
double sumC,sumP;		//sum of single call and put prices
ofstream output;

Random rnd;
rnd.SetRandom("Primes","seed.in");//initialization random_gen from files

//Point 1.
cout << endl << "Computing prices directly ..." << endl;
for(int i=0; i<N; i++){
    sumC=0.,sumP=0.;
    for(int j=0; j<L; j++){
	WT=rnd.Gauss(0.,1.);
	ST=S0*exp((r-sigma*sigma/2.)*T+sigma*WT*sqrt(T));
	if (ST > K)	{sumC+=exp(-r*T)*(ST-K);}	//in this case P=0
	else		{sumP+=exp(-r*T)*(K-ST);}	//(and sumC+=0.)
    }
    C.DefComp(i,sumC/L);
    P.DefComp(i,sumP/L);
}
C.StatErrProg(resC,erC);	//computation of res ed err, blocking method
P.StatErrProg(resP,erP);

//Point 2.
cout << "Computing prices, discretized GMB ..." << endl << endl;
C.SetUsed(0);P.SetUsed(0);
for(int i=0; i<N; i++){
    sumC=0.,sumP=0.;
    for(int j=0; j<L; j++){
        ST=S0;
	for(int k=0; k<Nt; k++){
	    WT=rnd.Gauss(0.,1.);		//sub-interval Wiener process
	    ST*=exp((r-sigma*sigma/2.)*(T/Nt)+sigma*WT*sqrt(T/Nt));
	}
	if (ST > K)	{sumC+=exp(-r*T)*(ST-K);}	//in this case P=0
	else		{sumP+=exp(-r*T)*(K-ST);}	//(and sumC+=0.)
    }
    C.DefComp(i,sumC/L);
    P.DefComp(i,sumP/L);
}
C.StatErrProg(resC_d,erC_d);	//computation of res ed err, blocking method
P.StatErrProg(resP_d,erP_d);

//writing results
output.open("res01.out");
output << N << endl << endl;
for (int i=0; i<N; i++)	{output << resC.GetComp(i) << "," << erC.GetComp(i) << "," << resP.GetComp(i) << "," << erP.GetComp(i) << "," << resC_d.GetComp(i) << "," << erC_d.GetComp(i) << "," << resP_d.GetComp(i) << "," << erP_d.GetComp(i) << endl;}
output.close();

return 0;
}
