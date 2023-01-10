#include "deps.h"
#include <stdio.h>
#include <bits/stdc++.h>

#include "globalconsts.h"
#include "functions_lib.h"

using namespace std;

const double eulergamma=0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495;

komplex::komplex()
{

}

komplex komplex::operator + (komplex rhs)
{
    komplex tmp;
    tmp.re=re+rhs.re;
    tmp.im=im+rhs.im;
    return tmp;
}

komplex komplex::operator * (komplex rhs)
{
    komplex tmp;
    tmp.re=re*rhs.re-im*rhs.im;
    tmp.im=im*rhs.re+re*rhs.im;
    return tmp;
}

komplex komplex::conj(komplex rhs)
{
    komplex tmp;
    tmp.re=rhs.re;
    tmp.im=-rhs.im;
    return tmp;
}

komplex:: ~komplex()
{

}

int factorial(int nn)
{
//    cout<<nn<<endl;
    if(nn==1 or nn==0)
    {
//        cout<<"i2"<<endl;
        return 1;
    }
 /*   if(nn>1)
    {
  //      cout<<"i3"<<endl;
        return nn*factorial(nn-1);
    }*/
    else
    {
        double s=1;
        for(int i=1;i<=nn;i++)
        {
            s=s*i;
        }
        return s;
    }
}

double log_factorial(double nn)
{
//    cout<<nn<<endl;
    if(nn==1. or nn==0.)
    {
//        cout<<"i2"<<endl;
        return 0;
    }
 /*   if(nn>1)
    {
  //      cout<<"i3"<<endl;
        return nn*factorial(nn-1);
    }*/
    else
    {
        double s=0;
        for(double i=1.;i<=nn;i++)
        {
            s=s+log(i);
        }
        return s;
    }
}

unordered_map<int,double> lfac;

double log_factorial2(double nn)
{
    int imin=0,imax=-1;
    if(lfac.size()<=nn)
    {
        imax=max(static_cast<int>(nn),1000);
        imin=max(0,static_cast<int>(lfac.size()));
    }

    for(int i=imin;i<imax;i++)
    {
        lfac[i]=log_factorial(static_cast<double>(i+1));
    }
    return lfac[static_cast<int>(nn-1)];
}

double harmonic_number(double k)
{
    double sum;
    for(double i=1;i<=k;i++)
    {
        sum=sum+ 1./i;
    }
    return sum;
}

double BesselJ_Re(double a, double r, double phi, int nmax)
{
    int A=static_cast<int>(a);
    double s=0;
    for(double m=0.;m<nmax;m++)
    {
  //      cout<<m<<" "<<A<<" "<<factorial(m)<<" "<<factorial(m+A)<<endl;
 //       s=s+pow(-1.,m)*pow(r/2.,2.*m+A)*cos((2.*m+A)*phi)/static_cast<double>(factorial(m)*factorial(m+A));
 s=s+pow(-1.,m)*exp((2.0*m+a)*log(r/2.0)-log_factorial2(m)-log_factorial2(m+a))*cos((2.*m+a)*phi);
    }
 //   cout<<"re"<<s<<endl;
    return s;
}

double BesselJ_Im(double a, double r, double phi, int nmax)
{
    int A=static_cast<int>(a);
    double s=0;
    for(int m=0;m<nmax;m++)
    {
 //       s=s+pow(-1.,m)*pow(r/2.,2.*m+A)*sin((2.*m+A)*phi)/static_cast<double>(factorial(m)*factorial(m+A));
  s=s+pow(-1.,m)*exp((2.0*m+a)*log(r/2.0)-log_factorial2(m)-log_factorial2(m+a))*sin((2.*m+a)*phi);

    }
   // cout<<"im"<<s<<endl;
    return s;
}

double BesselY_Re(double a, double r, double phi, int nmax)
{
    double s=(2.0/pi)*(eulergamma+log(r/2.0))*BesselJ_Re(a,r,phi,nmax)-(2.0/pi)*phi*BesselJ_Im(a,r,phi,nmax);
    int A=static_cast<int>(a);
    for(double k =0; k<=A-1;k++)
    {
        s=s-(1./pi)*(exp(log_factorial2(a-k-1)-log_factorial2(k)))*pow(r/2.0,2.0*k-a)*cos(phi*(2.0*k-a));
    }

    for(double k =0; k<=nmax;k++)
    {
        s=s-(1./pi)*pow(-1,k)*( (harmonic_number(k)+harmonic_number(k+a))*exp(-log_factorial2(k)-log_factorial2(k+a)))*pow(r/2.0,2.*k+a)*cos(phi*(2.0*k+a));
    }
    return s;
}

double BesselY_Im(double a, double r, double phi, int nmax)
{
    double s=(2.0/pi)*(eulergamma+log(r/2.0))*BesselJ_Im(a,r,phi,nmax)-(2.0/pi)*phi*BesselJ_Re(a,r,phi,nmax);
    int A=static_cast<int>(a);
    for(double k =0; k<=A-1;k++)
    {
        s=s-(1./pi)*(exp(log_factorial2(a-k-1)-log_factorial2(k)))*pow(r/2.0,2.0*k-a)*sin(phi*(2.0*k-a));
    }

    for(double k =0; k<=nmax;k++)
    {
        s=s-(1./pi)*pow(-1,k)*( (harmonic_number(k)+harmonic_number(k+a))*exp(-log_factorial2(k)-log_factorial2(k+a)))*pow(r/2.0,2.*k+a)*sin(phi*(2.0*k+a));
    }
    return s;
}

komplex BesselJ(double a , double r, double phi, int nmax)
{
    komplex tmp;
    tmp.re=BesselJ_Re( a ,  r,  phi,  nmax);
    tmp.im=BesselJ_Im( a ,  r,  phi,  nmax);
    return tmp;
}

komplex BesselY(double a , double r, double phi, int nmax)
{
    komplex tmp;
    tmp.re=BesselY_Re( a ,  r,  phi,  nmax);
    tmp.im=BesselY_Im( a ,  r,  phi,  nmax);
    return tmp;
}