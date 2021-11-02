#include"deps.h"

#include "TMDICE_lib.h"




int n=10000;
double gaussx[10000],gaussw[10000];
double N_F=3.0;
double N_c=3.0;

void fillgauss(double (&gaussx)[10000],double (&gaussw)[10000])
{
	gaussx[0]=0;
	gaussw[0]=0;
	for(int i=1;i<=n;i++)
	{
		gaussx[i]=cos(pi*(2.0f*i-1.0f)/(2.0f*n));
		gaussw[i]=sin(pi*(2.0f*i-1.0f)/(2.0f*n));
	}
}

double rand01(void)
{
return static_cast<  double>(rand()%RAND_MAX)/static_cast<  double>(RAND_MAX);
}

void fill_partitionfunction( fkt f, double vmin,double vmax,double vanz,map<double,double> &Finv, double &Fmax)
{
	double w=0;
	double dv=(vmax-vmin)/vanz;
	for(double v=vmin;v<=vmax-dv;v=v+dv)
	{
        double xx;
        w=integral(f,vmin,v,xx);
        Finv[w]=v;
	}
	Fmax=w;
}

void fill_partitionfunction_v2( fkt f, double vmin,double vmax,double vanz,map<double,double> &Finv, double &Fmax)
{

    double w=0;
    double dv=(vmax-vmin)/vanz;
    double fkt=0.001;

       for(double v=vmin;v<=vmax;v=v+dv)
    {
        double xx;
        w=w+integral(f,max(vmin,v-dv),v,xx);
        Finv[w]=v;
    }
    Fmax=w;
}

double select(fkt Finv, double Fmax,double wmin,double wmax)
{
	double r=Fmax*rand01();
	double w=Finv(r);
	w=max(wmin,w);
	w=min(wmax,w);
	return w;
}

double T(double x, int n)
{
	return cos(static_cast<double>(n)*acos(x));
}


