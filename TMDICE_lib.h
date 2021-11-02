#include"deps.h"

#ifndef TMDICE_lib	
#define TMDICE_lib

typedef double (*fkt)(double);

struct fvec
{
fkt q,g;
};

const long double pi =3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

	    const double C_F=4.0f/3.0f;
	    const double T_R=0.5f;
	    const double C_A=3.0f;


extern int n;
extern double gaussx[10000],gaussw[10000];
extern double N_F,N_c;

void fillgauss(double (&gaussx)[10000],double (&gaussw)[10000]);

template<typename T, typename ...Rest> double integral(  double f(  T x,Rest...rest), double a,  double b,T x, Rest...rest );

double rand01(void);

void fill_partitionfunction( fkt f, double vmin,double vmax,double vanz,map<double,double> &Finv, double &Fmax);

void fill_partitionfunction_v2( fkt f, double vmin,double vmax,double vanz,map<double,double> &Finv, double &Fmax);

double select(fkt Finv, double Fmax,double wmin,double wmax);


double T(double x, int n);



#include"TMDICE_lib.tpp"
#endif
