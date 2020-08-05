#include "Experiment.h"

#ifndef __Metropolis__
#define __Metropolis__

class Metropolis: public Experiment {

    protected:
    bool UseGauss;		//if 1 for Gaussian T, else uniform T (default)
    double HNow,HNew,Beta;	//AnnEnergy saved for efficiency, beta for Ann
    bool Ending;		//if we are in the SA final part
    
    //Internal functions
    void StdMove(DataVett,int ip);//effective single step of Markov Chain
    bool StdIfMove(DataVett,int ip);//verify if move accepted
    DataVett StdXnew(double) const;//Xnew from overall changing of X with T
    virtual double qRatio(const DataVett Xnew, int ip) const =0;//config new and all differ for particle ip (ip=-1 if total different config)

    virtual double AnnEnergy() const {return 0;};//energy eval for Sim Annealing
    virtual void PrintAcc() const {};//print some acceptance parameters
    double AnnqRatio() const;	//qRatio for Simulated Annealing

    public:
    //constructors
    Metropolis();
    Metropolis(Random* rnd);
    //destructor
    ~Metropolis(){}
    //Simulated Annealing (we have always to write specific AnnEnergy)
    void SimAnnealing(const DataVett beta, const DataVett bStep);

};

#endif
