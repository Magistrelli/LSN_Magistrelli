#include "Metropolis.h"

#ifndef __Hydrogen__
#define __Hydrogen__

class Hydrogen final: public Metropolis {

    private:
    int _n,_l,_m;		//state parameters
    const int iR=0,n_Props=1;	//radius index, number of observables
    double Param;		//T-matrix parameter

    //internal functions
    virtual void Move()	{StdMove(StdXnew(Param),-1);}
    virtual void Measure();
    virtual void Averages(int iblk, ofstream *OutRes);
    virtual double qRatio(const DataVett Xnew, int ip) const;
    
    public:
    //constructors
    Hydrogen(int nlm[3], const DataVett start, Random* rnd);
    Hydrogen(int[3],const DataVett,Random*, int usegauss);
    //destructor
    ~Hydrogen(){}
    //functions
    void SetParam(double param)	{Param=param;}

};


//spherical coordinates
DataVett Spheric(const DataVett XYZ);

#endif
