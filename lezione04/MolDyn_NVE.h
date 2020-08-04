#ifndef __MolDyn_h__
#define __MolDyn_h__

#include "Experiment.h"

class MolDyn: public Experiment {

    private:
    //thermodynamical state
    double Temp,Rho,Vol,Box,Rcut;
    //simulation
    double Delta;		
    //configuration (Pointer to TotD=DimSp*nPart-dim vectors)
    DataVett *Xold,*V;
    //parameters, observables
    const int iV=0,iK=1,iE=2,iT=3,iGr=4,nBins=100;
    const int n_Props=iGr+nBins;
    double BinSize;

    //Internal functions
    double Force(const int, const int) const;//Compute forces as -Grad_ip V(r)
    virtual void WriteInstant(int istep, ofstream* OutRes) const;
    void Move(DataVett*,double);	//Move particles with Verlet algorithm
    virtual void Move()	{Move(Xold,2.*Delta);}//stdMove is with Xold and 2dt
    virtual void Measure();			
    virtual void Averages(int iblk, ofstream* OutRes);
    double Pbc(const double r) const {return r-Box*rint(r/Box);}//side L=box
    DataVett Pbc(const DataVett&) const;
    void PrintInfo() const;		//Print simulation's info and params


    public:
    //constructor
    MolDyn();
    MolDyn(int nblk, Random* rnd);
    //distructor
    ~MolDyn();
    //functions
    int GetNBins() const	{return nBins;}
    double GetBinWidth() const	{return BinSize;}
    void Restart();			//restart simulation for warm-up
    void BoxScale()	{*X/=Box;}	//X in box scale
    void ConfFinalPlus() const;		//Write final configs, now and old

};

#endif
