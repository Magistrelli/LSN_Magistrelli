#include "Metropolis.h"

#ifndef __Variational_MC__
#define __Variational_MC__

#define hbar 1.		//Planck's costant

class PsiTrial final: public Metropolis { //here superposition of two gaussians (mu,sigma) and (-mu,sigma), but generalizable

    private:
    double Delta;		//T-matrix width
    Metropolis* Metro;		//for annealing of sampled function
    const unsigned int ParDim=2;//number of actual PsiT parameters (=X)
    const int iM=0,iS=1,n_Props=2;//obs are the parameters (mu,sigma)

    //internal functions
    virtual void Move();
    virtual void Measure();
    virtual void Averages(int iblk, ofstream *OutRes);
    virtual double qRatio(const DataVett,int) const	{return AnnqRatio();}
    virtual double AnnEnergy() const;
    virtual void PrintAcc() const;

    public:
    //constructors
    PsiTrial(Random* rnd);		//this if not use Metropolis for Psi
    PsiTrial(Random* rnd, double delta);//use this for SimAnneling
    //destructor
    ~PsiTrial(){}
    //functions
    void SetMetro(Metropolis* metro)	{Metro=metro;}
    double Eval(double x) const;	//function evaluation
    double Der2(double x) const;	//second derivative
    void SetPar(DataVett par)		{*X=par;}

};


class VarMC final: public Metropolis {

    private:
    double Mass=1.;		//Particle's mass
    PsiTrial* PsiT;		//trial wave function
    double Delta;		//T-matrix parameter
    //parameters, observables
    const int iE=0;		//energy index, number of observables
    const int iF=1,nBins=100;
    const int n_Props=iF+nBins;	//#obs, each bin of psi(x) is an obs
    double Xlim,BinSize;	//bin range: [-Xlim,Xlim]

    //internal functions
    double VExt(double x) const	{return pow(x,4)-5./2.*x*x;} //External V eval
    double Eloc(double x) const;//H*Psi/Psi, local energy
    virtual void Move()		{StdMove(StdXnew(Delta),-1);}
    virtual void Measure();
    virtual void Averages(int iblk, ofstream *OutRes);
    virtual double qRatio(const DataVett Xnew, int ip) const;
    
    public:
    //constructors
    VarMC(PsiTrial* psi, double xl, double del, Random* rnd);
    //destructor
    ~VarMC(){}
    //functions
    int GetNBins() const	{return nBins;}
    double GetBinWidth() const	{return BinSize;}

};

#endif
