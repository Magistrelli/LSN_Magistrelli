#include "Metropolis.h"
#include "TSPpath.h"

#ifndef __SimAnnTSP__
#define __SimAnnTSP__

//ATTENTION! here X from Metropolis and Gene from Path are the same thing!!

class MetroPath: public Metropolis, public Path {

    private:
    int StepMax;	//SA #step without emproved best path to end program
    int Istep,EndStep;	//progressive steps counters
    bool Square;	//if 1 cities in a square, else on a circumference
    string* Loc;	//usefull for output file
    double Box;		//square 2*side or circumference radius
    const int iL=0,n_Props=1;//obs is only path's lenght
    double MinLenght;

    //internal functions
    void SubMove();    		//single mutation trial
    virtual void Move();	//in every step try all the mutation operators
    virtual void Measure(){}	//useless, don't use blocking method (min search)
    virtual void Averages(int iblk, ofstream *OutRes){}//useless
    virtual double qRatio(const DataVett,int) const	{return AnnqRatio();}
    virtual double AnnEnergy() const			{return Loss;}
    virtual void PrintAcc() const;
    void WriteLength() const;	//write path's progressive lenght

    public:
    //constructors
    MetroPath(Random* rnd);
    //destructor
    ~MetroPath(){}

};

#endif
