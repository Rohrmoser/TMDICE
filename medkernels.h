#include"deps.h"
#include"TMDICE_lib.h"
#include"decor.h"
#include "globalconsts.h"
#include "functions_lib.h"

#include"TMDICE_base.h"
#include"BDIM.h"
#ifndef MEDKERNELS	
#define MEDKERNELS
extern map<double,map<double,fkt>>K;
extern map<double,map<double,fkt2>>Kt_stat,Kt_bj, Kt_exp,Kt;
extern map<double,map<double,fkt>>K_u;

extern map<double,map<double,fkt>>f;

extern map<double,map<double,map<double,double>>> phi_u_inv;
extern map<double,map<double,double>> phi_u_max;

double WW0(double R);
double WW1(double R);
double WW2(double R);
double WW3(double R);
#endif
