#include "random.h"
#include "Vettore.h"
#include <iomanip>

#ifndef __TSPpath__
#define __TSPpath__

#define Nspace 2	//cities on a plane

class Path {		//base single chromosome stuff

    protected:
    int Ncit;		//#cities to visit (in Gene will be 0,1,..,Ncit-1)
    DataVett *Pos;	//cities positions
    Random* _rnd;
    int* Gene;		//pointer to chromosome's first gene (GA's view)
    double Loss;	//path's cost
    int PathIndex;	//index of path in a possible group of them
    const int NMut=4;	//number of implemented mutation operators
    bool Flag;		//flag for correct evolution check

    //internal functions
    void MutSwitch();		//switch two random cities
    void MutShift();		//shift rnd# of contiguous cities at rnd dist
    void MutPermut();		//switch rnd# of contiguous cities at rnd dist
    void MutInvert();		//invert order of rnd# cities at rnd pos
    void SwitchGene(int,int);	//switch two genes' values
    double Distance(int,int) const;//distance between two cities
    void RndGen();		//first generation random generated
    void WritePath(ofstream&) const;//write actual path

    public:
    //constructors
    Path();
    Path(DataVett* pos, Random* rnd);
    //destructor
    ~Path() {delete[] Gene;}
    //values access
    int GetGene(int i) const	{return Gene[i];}
    void SetGene(int i, int a)	{Gene[i]=a;}
    double GetLoss() const	{return Loss;}
    void SetLoss(double loss)	{Loss=loss;}
    bool GetFlag() const	{return Flag;}
    int GetNcit() const		{return Ncit;}
    //functions
    void Check();		//check if every city appears 1! time
    void LossCost();		//loss function (total distance) eval
    void WriteStart(bool,double,string) const;//starting path in default file
    void WriteBest(string) const;//write actual path in best path default file

};


//external function
DataVett GenRndPos(int,bool,double,Random*);	//generate cities' random pos


#endif
