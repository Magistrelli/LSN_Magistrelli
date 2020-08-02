#ifndef __Vettore_h__
#define __Vettore_h__

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

class Vett {

	protected:
	unsigned int _N;
	unsigned int _used;
	double *_v;

	public:
	//costruttori
	Vett();
	Vett(unsigned int N);
	Vett(const Vett& vett);
	//distruttore
	~Vett();
	//ridefinizione operatori
	Vett& operator=(const Vett& vett);
	//definisce PER LA PRIMA VOLTA la componente k-esima
	void DefComp(unsigned int k, double value);
	//modifica la componente k-esima
	void SetComp(unsigned int k, double value);
	//set _used
	void SetUsed(unsigned int used)			{_used=used;}
	//accesso ai valori
	unsigned int GetDim() const			{return _N;}
	unsigned int GetUsed() const			{return _used;}
	double GetComp(unsigned int k) const;
	double* GetV() const				{return _v;}
	//stampa valori
	void Print() const;
	void PrintDebug() const;

};

class DataVett: public Vett {

	public:
	//costruttori
	DataVett();
	DataVett(unsigned int N);
	DataVett(unsigned int N,ifstream&);
	DataVett(const DataVett& vett);
	//distruttore
	~DataVett();
	//ridefinizione operatori
	DataVett& operator+=(const DataVett& vett);
	DataVett& operator*=(const double a);
	DataVett& operator/=(const double a);
	//funzioni
	void SumComp(unsigned int k, double value);//v[k]+=value
	void Read(ifstream&);			//inserisci dati da file
	void Print(ofstream&) const;		//stampa valori su file
	//res&err with blocking method (progressive)
	void StatErrProg(DataVett& res, DataVett& err) const;
};

//ridefinizione operatori
DataVett operator*(const double a, const DataVett& vett);	//prodotto per uno scalare
DataVett operator/(const DataVett& vett, const double a);
DataVett operator*(const DataVett& v1, const DataVett& v2);	//prodotto scalare
DataVett operator+(const DataVett& v1, const DataVett& v2);	//somma di vettori
DataVett operator-(const DataVett& v1, const DataVett& v2);

#endif
