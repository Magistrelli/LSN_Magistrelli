#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

#ifndef __Random__
#define __Random__

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // initializator
  void SetRandom(int* , int, int);
  void SetRandom(string, string);
  void SaveSeed();
  //gen from uniform distributions
  double Rannyu(void);
  double Rannyu(double min, double max);
  //gen from various distributions
  double Gauss(double mean, double sigma);
  double Exp(double lambda);
  double Lorentz(double mu, double Gamma);
  //gen uniformly an angle in [-pi,pi]
  double Angle();
};

#endif // __Random__
