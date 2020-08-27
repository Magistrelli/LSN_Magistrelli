#include "Experiment.h"

#ifndef __Metropolis__
#define __Metropolis__

class Metropolis: public Experiment {

    protected:
    bool UseGauss;	//if 1 for Gaussian T, else uniform T (default)
    double HNow,HNew,Beta;//for Simulated Annealing, Energy saved for efficiency
    bool Ending;		//if we are in the SA final part
    
    //Internal functions (config new and now differ for particle ip; ip=-1 if total different config)
    void StdMove(DataVett,int ip);//standard single step of Markov Chain
    bool StdIfMove(DataVett,int ip);//verify if move accepted
    DataVett StdXnew(double) const;//Xnew from overall changing of X with T
    virtual double qRatio(const DataVett Xnew, int ip) const =0;

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
