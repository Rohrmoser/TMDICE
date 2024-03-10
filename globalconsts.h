#include "deps.h"
#include "decor.h"

#ifndef GLOBALCONSTS
#define GLOBALCONSTS

const long double pi =3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

extern double gevfm;

extern double xeps,xmin;
extern double anz;
extern double qmin;
extern double qmax;
extern double qanz;
extern double u_min;
extern double u_max;
extern double u_anz;


extern double colorfact,ca,cf;
extern double tanz;
extern int time_mode; //time_mode=0...infinite constant medium, 1...static medium of length tmax

extern double pmax;
extern int dinorm;

extern double tmin,tmax,taumin,taumax, qhat,qhatgev, emax, alphas, alphabar,nc, ndens, inittyp,W_scat, md,tstar,tstargev,Temp,x0,kt0,typ0,R; 
extern double theta0orR,Q2max,xminvle,xepsvle, anzvle, Q2min,ao_mode, checkdla,vlemode,simplifythkt,setbias,rotate,p1x,p1y,p1z,bdimmode,vle2mode;
extern double scatflag1,ktsplitflag1,time_mode1,dlla;

extern int ktsplitflag, scatflag;

extern std::map<std::string,double> invars;


const int M=300;
const double m=static_cast<double>(M);

const double ee=1*pow(10,0),ee_u=7*pow(10,-5);

extern double cqq_u[M], cqg_u[M],cgg_u[M];

const int tanzmax=1*pow(10,4);//Maximal number of grid points. 
extern double c_t[3][3][tanzmax][M],c1_t[3][3][tanzmax][M],c2_t[3][3][tanzmax][M],c3_t[3][3][tanzmax][M];

extern std::map<double,std::map<double,double>> typ2;

void readTMDICEparameters(std::string snom);

void readTMDICEparameters(std::map<std::string,double>iv1);

extern std::map<double,double> up_lim_phi_time;

#endif
