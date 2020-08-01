#include "Vettore.h"

//VETTORE
//costruttori
Vett::Vett(): _N(0), _used(0){
	_v=NULL;
}
Vett::Vett(unsigned int N){
	_N=N;
	_used=0;
	_v = new double[N];
	for (unsigned int i=0; i<N; i++) _v[i]=0.;
}
Vett::Vett(const Vett& vett){
	_N=vett._N;
	_used=vett._used;
	_v = new double[_N];
	for (unsigned int i=0; i<_N; i++) _v[i]=vett._v[i];
}
//distruttore
Vett::~Vett(){
	delete[] _v;
}

//ridefinizione operatori
Vett& Vett::operator=(const Vett& vett){
	_N=vett._N;
	_used=vett._used;
	if (_v) delete[] _v;
	_v = new double[_N];
	for (unsigned int i=0; i<_N; i++) _v[i]=vett._v[i];
	return *this;
}

//definisce PER LA PRIMA VOLTA la componente k-esima
void Vett::DefComp(unsigned int k, double value){
	if (k<_N) {
		_v[k]=value;
		_used++;
	}
	else {
		cerr << endl << "ATTENZIONE!! DefComp: indice fuori dai limiti!!!" << endl;
		exit(-1);
	}
}
//modifica la componente k-esima
void Vett::SetComp(unsigned int k, double value){
	if (k<_used) _v[k]=value;
	else if (k<_N)
		cout << endl << "ATTENZIONE!! SetComp: indice superiore a used!!!" << endl;
	else {
		cerr << endl << "ATTENZIONE!! SetComp: indice fuori dai limiti!!!" << endl;
		exit(-1);
	}
}


//restituisce il valore della componente k-esima
double Vett::GetComp(unsigned int k) const{
	if (k<_used) return _v[k];
	else if (k<_N)
		cout << endl << "ATTENZIONE!! GetComp: indice superiore a used!!!" << endl;
	else {
		cerr << endl << "ATTENZIONE!! GetComp: indice fuori dai limiti!!!" << endl;
		exit (-1);
	}
return _v[k];
}
//restituisce il puntatore a al vettore
double* Vett::GetV() const{
	double *vett = new double [_used];
	for (unsigned int i=0; i<_used; i++)	{vett[i]=GetComp(i);}
	return vett;
}


//stampa valori
void Vett::Print() const{
	for (unsigned int i=0; i<_N; i++)
		cout << i << ")\t" << _v[i] << endl;
}
void Vett::PrintDebug() const{
	cout << "dimensione del vettore: " << _N << endl;
	cout << "used: " << _used << endl;
	cout << "componenti del vettore:" << endl;
	for (unsigned int i=0; i<_N; i++)
		cout << i << ")\t" << _v[i] << endl;
}



//VETTOREDATI
//costruttori
DataVett::DataVett(): Vett(){};
DataVett::DataVett(unsigned int N): Vett(N){};
DataVett::DataVett(unsigned int N, ifstream& input): Vett(N){
	Read(input);
}
DataVett::DataVett(const DataVett& vett): Vett(vett){};
//distruttore
DataVett::~DataVett(){}

//inserisci dati da file
void DataVett::Read(ifstream& input){
    if (input.fail())
	cout << endl << "ATTENZIONE!!! Problema apertura file!!!" << endl;
    else{
	while(!input.eof() && _used<_N){
		input>>_v[_used];
		_used++;
	}
//	if (!(*input).eof())	cout << endl << "ATTENZIONE! Il file non è stato letto completamente" << endl;
    }
}



//risultato con incertezza blocking method (progressive)
//risultato come media deelle singole misure, incertezza come stdDev della media
void DataVett::StatErrProg(DataVett& res, DataVett& err) const{
    double sum=0.,sum2=0.,ave,av2,sigma;
    for(unsigned int i=0; i<GetDim(); i++){
	sum+=GetComp(i);
	sum2+=pow(GetComp(i),2);
	ave=sum/(i+1.);
	av2=sum2/(i+1.);
	sigma=sqrt((i+1.)/i*(av2-pow(ave,2)));
	res.DefComp(i,ave);
	if(i==0){err.DefComp(i,0.);}
	else	{err.DefComp(i,sigma/sqrt(i));}
    }
}

//stampa valori su file
void DataVett::Print(ofstream& output) const{
	for (unsigned int i=0; i<_used; i++)
		output << i << ")\t" << _v[i] << endl;
}
