#include "GeneticTSP.h"

int main (int argc, char *argv[]){

//Variables
Random* rnd = new Random();
rnd->SetRandom("Primes","seed.in");	//initialization random_gen from files
Population Pop(rnd);
//int ncit=Pop.GetNcit();
//Chromo *Chr;

//simulation
/*//with this code you can test any of the mutation singularly put its prob = 1 and every other prob = 0, Ncit=8 and Nchr=5
for(int i=0; i<Pop.GetNchr(); ++i){
    Chr=Pop.GetChromo(i);
    cout<<endl<< i+1 << "-th Chromosome:"<<endl;
    for(int j=0; j<ncit; ++j){
	cout<<"  "<< Chr->GetGene(j)+1;
    }
    cout<<endl;
}
cout<<endl<<endl<<"After"<<endl;
for(int i=0; i<Pop.GetNchr(); ++i){
    Chr=Pop.GetChromo(i);
    Chr->Mutation();
    cout<<endl<< i+1 << "-th Chromosome:"<<endl;
    for(int j=0; j<ncit; ++j){
	cout<<"  "<< Chr->GetGene(j)+1;
    }
    cout<<endl;
}*/

Pop.Evolution();		//comment if mutation check
cout<<endl;

return 0;
}
