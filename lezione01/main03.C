#include "Vettore.h"
#include "random.h"


int main (int argc, char *argv[]){

//Variables
double L=0.5,D=1.;		//needle lenght and distance between lines
int Nlines=int(1e4);		//number of lines
int M=int(1e7),N=int(2e2);	//total number of throws and blocks
int Mb=int(M/N);		//total number of throws in each block
int count,hit;			//progressive number of throws and success
double mean,angle;		//random mean point and inclination of the needle
double ymin,ymax;		//vertical projection of needle's extremes
int which;			//line identifier
DataVett Pi(N),res(N),err(N);	//single measure, progressive averages and errors
ofstream output;

Random rnd;
rnd.SetRandom("Primes","seed.in");//initialization random_gen from files

//Computing
for(int i=0; i<N; i++){
    count=0;
    hit=0;
    while(count<Mb){
	mean=rnd.Rannyu(0.,D*(Nlines-1.));
	angle=fabs(rnd.Angle());//positive or not give the same physical situation
	count++;
	ymax=mean+(L/2.)*sin(angle);
	ymin=mean-(L/2.)*sin(angle);
	which=int(ymin/D)+1;
	if(int(ymax/D)==which)	{hit++;}
    }
    Pi.DefComp(i,(2.*L*Mb)/(hit*D));
}
Pi.StatErrProg(res,err);	//computation of res ed err, blocking method

output.open("res03.out");
output << N << endl << endl;
for (int i=0; i<N; i++)	{output<< res.GetComp(i) << "," << err.GetComp(i) <<endl;}
output.close();
return 0;
}
