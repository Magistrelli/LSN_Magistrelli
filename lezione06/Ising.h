#include "Metropolis.h"

#ifndef __ISING__
#define __ISING__

class Ising: public Metropolis {

    private:
    //thermodynamical state
    double Temp,Beta,J,h;
    //simulation
    bool UseMetro;				//1 for true, 0 for false
    //parameters, observables
    const int iU=0,iC=1,iM=2,iX=3,n_Props=4;	//Energy, Heat capacity, Magnetization, Magnetic susceptibility index, Number of observables

    //Internal functions
    virtual void Move();
    virtual double qRatio(const DataVett Snew, int ip) const;
    void Gibbs(int o);
    virtual void Measure();
    virtual void Averages(int iblk, ofstream *OutRes);
    int Pbc(int) const;

    public:
    //constructors
    Ising(Random*);
    //destructor
    ~Ising(){}
    //functions
    bool GetMetro() const	{return UseMetro;}
    double Geth() const		{return h;}
    void SetTemp(double T_in)	{Temp=T_in,Beta=1./T_in;}
    void PrintInfo() const;

};

#endif
