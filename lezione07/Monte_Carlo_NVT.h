#include "Metropolis.h"

#ifndef __NVT__
#define __NVT__

class CanonicEns: public Metropolis {

    private:
    //thermodynamical state
    double Temp,Beta,Rho,Vol,Box,Rcut;
    //simulation
    double Delta;			//Uniform-T width
    //parameters, observables
    const int iV=0,iW=1;		//Potential energy, virial
    double Vtail,Ptail;			//Tail corrections
    //Measurement of g(r)
    const int iGr=2,nBins=100;
    const int n_Props=iGr+nBins;	//#obs, each bin of g(r) is an obs
    double BinSize;

    //Internal functions
    virtual void WriteInstant(int istep, ofstream* OutRes) const;
    virtual void Move();
    virtual void Measure();
    virtual void Averages(int iblk, ofstream *OutRes);
    virtual double qRatio(const DataVett Xnew, int ip) const;
    double Pbc(const double r) const {return r-Box*rint(r/Box);}//side L=box
    void PrintInfo() const;
    //Observables' values per particle
    double Epot(const double v) const	{return v/double(nPart)+Vtail;}
    double Virial(const double w) const	{return w/double(nPart)+Ptail;}
    double Pressure(const double w) const{return Rho*Temp+(w+double(nPart)*Ptail)/Vol;}

    public:
    //constructors
    CanonicEns(Random*);
    //destructor
    ~CanonicEns(){}
    //functions
    int GetNBins() const	{return nBins;}
    double GetBinWidth() const	{return BinSize;}
    void BoxScale()		{*X/=Box;}	//X in box units

};

#endif
