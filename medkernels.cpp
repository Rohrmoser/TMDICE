#include"deps.h"
#include"TMDICE_lib.h"
#include"decor.h"
#include "globalconsts.h"
#include "functions_lib.h"

#include"TMDICE_base.h"
#include"BDIM.h"
using namespace std;

double kqq(double z){return 1*0.5*cf*(1.0+z*z)*pow(1.0-z,-1)*sqrt(z*ca+pow(1.0-z,2)*cf)/sqrt(z*(1.0-z));}
double kqg(double z){double wert=1*N_F*T_R*(z*z+pow(1.0-z,2)) *sqrt(cf-z*(1.0-z)*ca)/sqrt(z*(1.0-z));return wert/(2.);}
double kgq(double z){return 1*0.5*cf*(1.0+pow(1.0-z,2))*(1.0/z)*sqrt((1.0-z)*ca+z*z*cf)/sqrt(z*(1.0-z));}
double kgg(double z){double wert=pow(ca,1.5)*pow(1.0-z+z*z,2.5)*pow(z*(1.0-z),-1.5); return wert/2.;}
map<double,map<double,fkt>>K={{1.0, {{1.0,kqq},{2.0,kgq}}},{-1.0, {{-1.0,kqq},{2.0,kgq}}},{2.0, {{1.0,kqg},{2.0,kgg}}}};


double kqq_stat(double z,double t){double kz=sqrt((z*ca+(1.0-z)*(1.0-z)*cf)/(z*(1-z))),tau=max(0.,(t)/(/*alphabar*/tstar) ); 
return 1.000000*kqq(z)*max(0.,-sin(kz*tau)/cosh(kz*tau)+tanh(kz*tau))/(1.+cos(kz*tau)/cosh(kz*tau));}
double kqg_stat(double z, double t){double kz=sqrt((cf-z*(1.-z)*ca)/(z*(1-z))),tau=max(0.,(t)/(/*alphabar*/tstar) ); 
return 1.000000*kqg(z)*max(0.,-sin(kz*tau)/cosh(kz*tau)+tanh(kz*tau))/(1.+cos(kz*tau)/cosh(kz*tau));}
double kgq_stat(double z,double t){double kz=sqrt(((1.-z)*ca+z*z*cf)/(z*(1-z))),tau=max(0.,(t)/(/*alphabar*/tstar) ); 
return 1.000000*kgq(z)*max(0.,-sin(kz*tau)/cosh(kz*tau)+tanh(kz*tau))/(1.+cos(kz*tau)/cosh(kz*tau));}
double kgg_stat(double z, double t){double kz=sqrt(ca*(1-z+z*z)/(z*(1-z))),tau=max(0.,(/*ca*/t/*-tmin*/)/(/*alphabar*/tstar) ); 
 return /*pow(ca,-1.5)*/1*kgg(z)*max(0.,-sin(kz*tau)/cosh(kz*tau)+tanh(kz*tau))/(1.+cos(kz*tau)/cosh(kz*tau));}
map<double,map<double,fkt2>>Kt_stat={{1.0, {{1.0,kqq_stat},{2.0,kgq_stat}}},{-1.0, {{-1.0,kqq_stat},{2.0,kgq_stat}}},{2.0, {{1.0,kqg_stat},{2.0,kgg_stat}}}};

double kqq_exp(double z,double t){double kz=sqrt((z*ca+(1.0-z)*(1.0-z)*cf)/(z*(1-z))),tau=max(0.,(t)/(alphabar*tstar) ); return 1.000000*kqq(z)*max(0.,-sin(kz*tau)/cosh(kz*tau)+tanh(kz*tau))/(1.+cos(kz*tau)/cosh(kz*tau));}
double kqg_exp(double z, double t){double kz=sqrt((cf-z*(1.-z)*ca)/(z*(1-z))),tau=max(0.,(t)/(alphabar*tstar) ); return 1.000000*kqg(z)*max(0.,-sin(kz*tau)/cosh(kz*tau)+tanh(kz*tau))/(1.+cos(kz*tau)/cosh(kz*tau));}
double kgq_exp(double z,double t){double kz=sqrt(((1.-z)*ca+z*z*cf)/(z*(1-z))),tau=max(0.,(t)/(alphabar*tstar) ); return 1.000000*kgq(z)*max(0.,-sin(kz*tau)/cosh(kz*tau)+tanh(kz*tau))/(1.+cos(kz*tau)/cosh(kz*tau));}
double kgg_exp(double z, double t){double kz=sqrt(/*ca*/(1-z+z*z)/(z*(1-z))),tau=max(0.,(ca*t/*-tmin*/)/(/*alphabar*/tstar) ); 
double r=kz*tau*sqrt(2.),xi=-pi/4.0;
int nmax=100;
double rj0=BesselJ_Re(0,r,xi,nmax),ij0=BesselJ_Im(0,r,xi,nmax),rj1=BesselJ_Re(1,r,xi,nmax),ij1=BesselJ_Im(1,r,xi,nmax);
double timepart=(-rj1*rj0-ij0*ij1+rj1*ij0-ij1*rj0)/(rj0*rj0+ij0*ij0);
return 1.0*pow(ca,-0.5)*kgg(z)*timepart;}
map<double,map<double,fkt2>>Kt_exp={{1.0, {{1.0,kqq_exp},{2.0,kgq_exp}}},{-1.0, {{-1.0,kqq_exp},{2.0,kgq_exp}}},{2.0, {{1.0,kqg_exp},{2.0,kgg_exp}}}};

double kqq_bj(double z,double t){double kz=sqrt((z*ca+(1.0-z)*(1.0-z)*cf)/(z*(1-z))),tau=max(0.,(t)/(alphabar*tstar) ); return 1.000000*kqq(z)*max(0.,-sin(kz*tau)/cosh(kz*tau)+tanh(kz*tau))/(1.+cos(kz*tau)/cosh(kz*tau));}
double kqg_bj(double z, double t){double kz=sqrt((cf-z*(1.-z)*ca)/(z*(1-z))),tau=max(0.,(t)/(alphabar*tstar) ); return 1.000000*kqg(z)*max(0.,-sin(kz*tau)/cosh(kz*tau)+tanh(kz*tau))/(1.+cos(kz*tau)/cosh(kz*tau));}
double kgq_bj(double z,double t){double kz=sqrt(((1.-z)*ca+z*z*cf)/(z*(1-z))),tau=max(0.,(t)/(alphabar*tstar) ); return 1.000000*kgq(z)*max(0.,-sin(kz*tau)/cosh(kz*tau)+tanh(kz*tau))/(1.+cos(kz*tau)/cosh(kz*tau));}
double kgg_bj(double z, double t){double kz=sqrt(/*ca*/(1-z+z*z)/(z*(1-z))),tau=max(0.,(ca*t/*-tmin*/)/(/*alphabar*/tstar) ); 
double r=kz*tau*sqrt(2.),xi=-pi/4.0;
int nmax=10;
double rj0=BesselJ_Re(0,r,xi,nmax),ij0=BesselJ_Im(0,r,xi,nmax),rj1=BesselJ_Re(1,r,xi,nmax),ij1=BesselJ_Im(1,r,xi,nmax);
double timepart;
return 1.0*pow(ca,-0.5)*kgg(z)*timepart;}
map<double,map<double,fkt2>>Kt_bj={{1.0, {{1.0,kqq_bj},{2.0,kgq_bj}}},{-1.0, {{-1.0,kqq_bj},{2.0,kgq_bj}}},{2.0, {{1.0,kqg_bj},{2.0,kgg_bj}}}};

map<double,map<double,fkt2>>Kt;

double k_u_qq(double z){return sin(z)*exp(-z);;}
double k_u_qg(double z){return sin(z)*exp(-z);;}
double k_u_gq(double z){return sin(z)*exp(-z);;}
double k_u_gg(double z){return sin(z)*exp(-z);}
map<double,map<double,fkt>>K_u={{1.0, {{1.0,k_u_qq},{2.0,k_u_gq}}},{-1.0, {{-1.0,k_u_qq},{2.0,k_u_gq}}},{2.0, {{1.0,k_u_qg},{2.0,k_u_gg}}}};

double fgg(double z){return 1*((1.-z)+z*z)*(ca);}
double fqg(double z){return cf-z*(1.-z)*ca;}
double fgq(double z){return (1.-z)*ca+z*z*cf;}
double fqq(double z){return z*ca+(1.-z)*(1.-z)*cf;}
map<double,map<double,fkt>>f={{1.0, {{1.0,fqq},{2.0,fgq}}},{-1.0, {{-1.0,fqq},{2.0,fgq}}},{2.0, {{1.0,fqg},{2.0,fgg}}}};

map<double,map<double,map<double,double>>> phi_u_inv;
map<double,map<double,double>> phi_u_max;

double WW0(double R){return 0;}
double WW1(double R){return qmin*qmin/(1.0-R);}
double WW2(double R){return md*md/(pow(1.0+md*md/(qmin*qmin),1.0-R)-1.0);}
double WW3(double R){return md*md/(pow(0.5*(1.0+md*md/(qmin*qmin) ),1.0-R)-1.0);}