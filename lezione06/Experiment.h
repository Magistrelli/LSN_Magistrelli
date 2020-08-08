#include "random.h"
#include "Vettore.h"
#include <iomanip>

#ifndef __Experiment__
#define __Experiment__

class Experiment {

    protected:
    //simulation
    int nPart,DimSp,TotD;//number of degrees of freedom (#DoF) and space's dim
    int nBlk,nStep;	//number of blocks and steps for each block
    int nBurn;		//number of burn-in steps
    int accepted,tries;	//accumulator for some acceptance estimation
    double TotA;	//sum of all blocks' acceptance estimator
    Random* Rnd;
    int iPrint=int(1e7);//each iPrint steps print #steps done
    bool _restart;	//1 for true, 0 for false; deafaul: false
    //configuration (Pointer to TotD=DimSp*nPart-dim vectors)
    DataVett *X;	//will have dimension = #DoF
    //observables
    int m_Props;	//will be number of obs with progressive res/err graph
    DataVett *Walker;	//accumulators for each observables
    DataVett *GlobAv,*GlobAv2,*Err;//blocking method

    //internal functions
    void InitVett();		//dynamic allocation of vectors
    void Reset(int iblk);	//Reset block accumulators
    void DoAverages(DataVett stima, int iblk) const;
    void WriteAverages(DataVett stima, int iblk, ofstream* OutRes) const;
    virtual void WriteInstant(int istep, ofstream* OutRes) const {};//for InstantValues, to be implemented in derived class
    void WriteConf(ofstream&,const int,const int,int) const;//last blocks conf
    void ConfXYZ(const int nconf) const;//Write configuration in .xyz format
    virtual void Move() =0;		//single specific algorithm move
    virtual void Measure() =0;		//Properties measurement
    virtual void Averages(int blk, ofstream *OutRes) =0;

    public:
    //constructors
    Experiment();
    Experiment(Random* rnd);
    //destructor
    ~Experiment();
    //values access
    int GetDim() const	{return nPart;}
    int GetL() const	{return nStep;}
    int GetBlk() const	{return nBlk;}
    //manual setting
    void SetBlkSt(int nblk, int nstep)	{nBlk=nblk,nStep=nstep;}
    void SetBurn(int nburn)		{nBurn=nburn;}
    void SetiPrint(int iprint)		{iPrint=iprint;}
    //functions
    double Acceptance() const	{return TotA/nBlk;}
    void Restart(const DataVett start);	//Re-initialize
    double* GetResErr(int iObs) const;	//final results and errors
    void ConfFinal() const;		//write final config (deafaul file)
    void ConfFinal(ofstream&) const;	//write final configuration
    //Simulation
    void DoSampling(ofstream* OutRes, ofstream& OutConf, int nnconf);
    void InstantValues(ofstream* OutRes, int nobs, int iout, int iconf);
    void BurnIn(ofstream& OutBurn, bool out);	//burn in, write Walker[0]
    void BurnIn();				//write in default file

};


//statistical uncertanties for Blocking Method
double Error(double,double,int);

#endif
