#include"deps.h"
#include"globalconsts.h"


#ifndef TMDICE_lib	
#define TMDICE_lib

typedef double (*fkt)(double);
typedef double (*fkt2)(double,double);


struct fvec
{
fkt q,g;
};


const double C_F=4.0f/3.0f;
const double T_R=0.5f;
const double C_A=3.0f;


extern int n;
extern double gaussx[10000],gaussw[10000];
extern double N_F,N_c;

void fillgauss(double (&gaussx)[10000],double (&gaussw)[10000]);

template<typename T, typename ...Rest> double integral(  double f(  T x,Rest...rest), double a,  double b,T x, Rest...rest );

double rand01(void);


template<typename T, typename ...Rest> void fill_partitionfunction_v2(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,map<double,double> &Finv, double &Fmax,T x, Rest...rest );

template<typename T, typename ...Rest> void fill_partitionfunction_v5(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,unordered_map<int,double> &Finv, double &Fmax,T x, Rest...rest );

template<typename T, typename ...Rest> void fill_partitionfunction_v6(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,unordered_map<int,double> &Finv, double &Fmax,T x, Rest...rest );

template<typename T, typename ...Rest> void fill_partitionfunction_v7(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,unordered_map<int,double> &Finv, double &Fmax,T x, Rest...rest );

template<typename T, typename ...Rest> void fill_partitionfunction_v8(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,unordered_map<int,double> &Finv, double &Fmax,T x, Rest...rest );

template<typename T, typename ...Rest> void fill_partitionfunction_v3(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,map<double,double> &Finv, double &Fmax,T x, Rest...rest );

template<typename T, typename ...Rest> void fill_partitionfunction_v4(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,map<double,double> &Finv, double &Fmax,T x, Rest...rest );


void fill_partitionfunction( fkt f, double vmin,double vmax,double vanz,map<double,double> &Finv, double &Fmax);

void fill_partitionfunction_v2( fkt f, double vmin,double vmax,double vanz,map<double,double> &Finv, double &Fmax);

void fill_partitionfunction_v3( fkt f, double vmin,double vmax,double vanz,map<double,double> &Finv, double &Fmax);


double select(fkt Finv, double Fmax,double wmin,double wmax);

double T(double x, int n);

vector<double> det_boundaries(fkt f, double y, double zmin,double zmax);

void make_c(map<double,double> f,double (&c)[M],double a, double b,double eee);

void make_3c(map<double,double>f,double (&c1)[M],double (&c2)[M],double (&c3)[M],double a,double b1,double b2,double b);

double approx_f(double x, double a, double b,double (&c)[M],double eee);

double approx_f(double x,double *c1,double *c2,double *c3,double a,double b1,double b2,double b,int M);

double interpolate_integer_inversegrid_1d(unordered_map<int,double> *phinvtmp,double xminus, double xplus, double r);

//double select(fkt Finv, double Fmax,double wmin,double wmax);

double select(fkt2 Finv, double Fmax,double wmin,double wmax,double t);




vector<double> det_boundaries(fkt f, double y, double zmin,double zmax);

vector<double> det_boundaries(fkt2 f, double y, double zmin,double zmax,double t);

void make_3c(map<double,double>f,double (&c1)[M],double (&c2)[M],double (&c3)[M],double a,double b1,double b2,double b);

//double approx_f(double x,double (&c1)[M],double (&c2)[M],double (&c3)[M],double a,double b1,double b2,double b);

void cheb_coeff_split(map<double,map<double,fkt2>>kernel,map<double,map<double,map<double,double>>>ph_max,map<double,map<double,map<double,map<double,double>>>>phi_inv,int f1, int f2
, double ymin,double (&c1)[3][3][tanzmax][M],double (&c2)[3][3][tanzmax][M],double (&c3)[3][3][tanzmax][M]
,map<double,map<double,map<double,double>>> (&b1),map<double,map<double,map<double,double>>> (&b2),map<double,map<double,map<double,double>>> (&b),double t);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include"TMDICE_lib.tpp"
#endif
