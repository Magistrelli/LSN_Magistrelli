#include "random.h"

Random :: Random(){}

Random :: ~Random(){}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

void Random::SetRandom(string fPrimes, string fSeed){
	int s[4],p1,p2;
	string property;

	ifstream inputPrimes(fPrimes);
	if (inputPrimes.is_open()){
		inputPrimes >> p1 >> p2;
	} else cerr << "PROBLEM: Unable to open " << fPrimes << endl;
	inputPrimes.close();
	
	ifstream inputSeed(fSeed);
	if (inputSeed.is_open()){
		while (!inputSeed.eof()){
			inputSeed >> property;
			if(property=="RANDOMSEED"){
				inputSeed >> s[0] >> s[1] >> s[2] >> s[3];
				SetRandom(s,p1,p2);
			}
		}
		inputSeed.close();
	} else cerr << "PROBLEM: Unable to open " << fSeed << endl;
	
return;
}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << "RANDOMSEED\t" << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Exp(double lambda) {
   double y=Rannyu();
   return -log(1.-y)/lambda;
}

double Random :: Lorentz(double mu, double Gamma) {
   double y=Rannyu();
   return mu+Gamma*tan(M_PI*(y-0.5));
}

double Random::Angle(){
   double x,y,hyp;
   do{	x=Rannyu(-1.,1.);
   	y=Rannyu(-1.,1.);
   	hyp=x*x+y*y;	}
   while (hyp>=1.);
   if(y<0.)	{return -acos(x/sqrt(hyp));}
   else		{return acos(x/sqrt(hyp));}
}
