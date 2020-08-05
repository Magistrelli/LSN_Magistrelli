#include "Hydrogen.h"

Hydrogen::Hydrogen(int nlm[3], const DataVett start, Random* rnd): Metropolis(rnd){
    _n=nlm[0],_l=nlm[1],_m=nlm[2];
    nPart=1,DimSp=3;			//1 particle, 3D space
    TotD=DimSp*nPart;
    m_Props=n_Props;
    Param=0.;

    InitVett();
    X->SetUsed(DimSp*nPart);
    Restart(start);			//initialization at start
}

Hydrogen::Hydrogen(int nlm[3], const DataVett start, Random* rnd , int usegauss): Hydrogen(nlm,start,rnd){
    UseGauss=usegauss;
}


//config new and now are totally different (ip=-1)
double Hydrogen::qRatio(const DataVett Xnew, int ip) const{
    if(ip!=-1)	{cerr<<endl<<endl<<"Attention! You should cange all the x,y,z particles' coordinates (ip should be -1)!"<<endl<<endl;}
    DataVett Sph1=Spheric(Xnew),Sph2=Spheric(*X);
    double r1=Sph1.GetComp(0),theta1=Sph1.GetComp(1),phi1=Sph1.GetComp(2);
    double r2=Sph2.GetComp(0),theta2=Sph2.GetComp(1),phi2=Sph2.GetComp(2);

    if ((_n==1) & (_l==0) & (_m==0)){return exp(-2.*(r1-r2));}
    else if ((_n==2) & (_l==1) & (_m==0)){return pow(r1*cos(theta1)/(r2*cos(theta2)),2)*exp(r2-r1);}
    else{
	cerr << endl << endl << "No rule to eval Hydro for n=" << _n << ", l=" << _l << " and m=" << _m << endl << endl;
	return 0;
    }
}


//Evaluate distance from the origin of the actual configuration
void Hydrogen::Measure(){
    double r2=0.;
    for(int i=0; i<TotD; i++)	{r2+=pow(X->GetComp(i),2);}
    Walker->SumComp(iR,sqrt(r2));
}
//Print results for current block
void Hydrogen::Averages(int iblk, ofstream *OutRes){
    DataVett stima(n_Props);
    stima.DefComp(iR,Walker->GetComp(iR)/double(nStep));
    DoAverages(stima,iblk);
    WriteAverages(stima,iblk,OutRes);
}





//EXTERNAL FUNCTIONS

//spherical coordinates
DataVett Spheric(const DataVett XYZ) {
    double x,y,z,r,theta=0.,phi=0.;
    DataVett Sph(3);
    int ir=0,iTh=1,iphi=2;
    if (XYZ.GetDim()!=3) {cerr<<endl<<endl<<"Hydrogen athom problem with DimSp = 3, given point of a "<< XYZ.GetDim() <<"-dimensional space!"<<endl<< endl;}
    x=XYZ.GetComp(0);
    y=XYZ.GetComp(1);
    z=XYZ.GetComp(2);
    
    r=sqrt(x*x+y*y+z*z);
    theta=acos(z/r);
    if(y<0.)	{phi=-acos(x/(r*sin(theta)));}
    else	{phi=acos(x/(r*sin(theta)));}

    Sph.DefComp(ir,r);
    Sph.DefComp(iTh,theta);
    Sph.DefComp(iphi,phi);
    return Sph;
}
