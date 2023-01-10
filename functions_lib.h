#include "deps.h"

#include "globalconsts.h"

using namespace std;

#ifndef FUNCTIONS_LIB
#define FUNCTIONS_LIB

class komplex
{
    public:
    double re,im;
    komplex();

    komplex operator + (komplex rhs);
    komplex operator * (komplex rhs);
    komplex conj(komplex rhs);
    virtual ~komplex();
};

//const double gamma=0.57721 56649 01532 86060 65120 90082 40243 10421 59335 93992 35988 05767 23488 48677 26777 66467 09369 47063 29174 67495;

int factorial(int n);

double log_factorial(double nn);

double log_factorial2(double nn);

double harmonic_number(int k);

double BesselJ_Re(double a, double r, double phi, int nmax=100);

double BesselJ_Im(double a, double r, double phi, int nmax=100);

double BesselY_Re(double a, double r, double phi, int nmax=100);

double BesselY_Im(double a, double r, double phi, int nmax=100);

komplex BesselJ(double a, double r, double phi, int nmax=100);

komplex BesselY(double a, double r, double phi, int nmax=100);

#endif