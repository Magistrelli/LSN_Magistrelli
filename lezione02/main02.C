#include "Vettore.h"
#include "random.h"


int main (int argc, char *argv[]){

//Variables
int M=int(5e4),N=int(1e2);	//number of total RWs and blocks
int L=int(M/N);			//number of RWs in each block
int Nst=int(1e2)+1;		//number of steps
double a=1.;			//lattice constant and steps' lenght
double x0=0.,y0=0.,z0=0.;	//starting point's coordinates
double x,y,z,r2;		//progressive point's position and dist^2
DataVett sum(Nst),Zero(Nst);	//sum all RW in each block, zero vect
DataVett Sum(Nst),Sum2(Nst);	//progressive sums for blocking method
DataVett res(Nst),err(Nst);	//means and errors for each time
DataVett resC(Nst),erC(Nst);	//continuum case
double sort,ave,ave2;		//random number, averages
double theta,phi;		//solid angle's variables
ofstream output;

Random rnd;
rnd.SetRandom("Primes","seed.in");//initialization random_gen from files
Zero.SetUsed(Nst);		//initialization of sums DataVett
Sum.SetUsed(Nst);
Sum2.SetUsed(Nst);

//Point 1.
for(int k=0; k<N; k++){		//repeat for each of N experiments
    sum=Zero;			//re-inizialization of sum for each exp
    for(int j=0; j<L; j++){	//for each experiment L random walks
	x=x0,y=y0,z=z0;
	r2=x*x+y*y+z*z;
	sum.SumComp(0,r2);
	for(int i=1; i<Nst; i++){	//creating a single RW
	    sort=rnd.Rannyu(0.,3.);	//from uniform distribution in [0,3] (equal probability for each direction of motion)
	    if(sort<1.){		//in [0,1] motion in x direction
		if (sort<0.5)	{x-=a;}	//negative motion
		else		{x+=a;}	//positive motion
	    } else if(sort<2.){		//in [1,2] motion in y direction
		if (sort<1.5)	{y-=a;}
		else		{y+=a;}
	    }				//in [2,3] motion in z direction
	    else if(sort<2.5)	{z-=a;}
	    else		{z+=a;}
	    r2=x*x+y*y+z*z;
	    sum.SumComp(i,r2);
	}
    } //ended a block, update global sums for each step
    for(int i=0; i<Nst; i++){
	ave=sum.GetComp(i)/L;		//sigle block's mean distance(i)^2
	Sum.SumComp(i,sqrt(ave));	//sum of |dist(i)| over blocks
	Sum2.SumComp(i,ave);		//sum for variance
    }
}
for(int i=0; i<Nst; i++){	//final results
    ave=Sum.GetComp(i)/N;	//mean value of |dist(i)|
    res.DefComp(i,ave);
    ave2=Sum2.GetComp(i)/N;	//mean value of dist(i)^2
    if (i==0){err.DefComp(0,0.);}
    else    {err.DefComp(i,sqrt((ave2-ave*ave)/(N-1.)));}//stdDev of the mean
}

//Point2.
Sum=Zero;			//re-initialization of sums
Sum2=Zero;
for(int k=0; k<N; k++){
    sum=Zero;
    for(int j=0; j<L; j++){
	x=x0,y=y0,z=z0;
	r2=x*x+y*y+z*z;
	sum.SumComp(0,r2);
	    for(int i=1; i<Nst; i++){
		rnd.SolidAngle(theta,phi);	//uniform sorted solid angle
		x+=a*sin(theta)*cos(phi);
		y+=a*sin(theta)*sin(phi);
		z+=a*cos(theta);
		r2=x*x+y*y+z*z;
		sum.SumComp(i,r2);
	}
    }
    for(int i=0; i<Nst; i++){
	ave=sum.GetComp(i)/L;
	Sum.SumComp(i,sqrt(ave));
	Sum2.SumComp(i,ave);
    }
}
for(int i=0; i<Nst; i++){
    ave=Sum.GetComp(i)/N;
    resC.DefComp(i,ave);
    ave2=Sum2.GetComp(i)/N;
    if (i==0){erC.DefComp(0,0.);}
    else     {erC.DefComp(i,sqrt((ave2-ave*ave)/(N-1.)));}
}

//writing results
output.open("res02.out");	//each line will represent a discrete time
output << Nst << "," << a << endl << endl;
for (int i=0; i<Nst; i++)	{output << res.GetComp(i) << "," << err.GetComp(i) << "," << resC.GetComp(i) << "," << erC.GetComp(i) << endl;}
output.close();
return 0;
}
