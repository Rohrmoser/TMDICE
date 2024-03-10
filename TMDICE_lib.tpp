#include<iostream>
#include<math.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include"deps.h"
#include"globalconsts.h"


#include "TMDICE_lib.h"
using namespace std;

	

   	template<typename T, typename ...Rest> double integral(  double f(  T x,Rest...rest), double a,  double b,T x, Rest...rest )
    	{
    	double iwert=0;
		if(b>a){
			for(int i=1;i<=n;i++)
			{
			iwert=iwert+
			f((a+b)/2.0f+((b-a)/2.0f)*gaussx[i], rest...)*gaussw[i];
			}
		}
	return iwert*pi/n*(b-a)/2.0f;
	}

	template<typename T, typename ...Rest> void fill_partitionfunction_v9(  double f(  T x,Rest...rest), double vmin,  
	double vmax, double vanz,unordered_map<int,double> &F,unordered_map<int,double> &Finv, double &Fmax,T x, Rest...rest )
	{
		double (*fpt)(T, Rest...)=(f);
		
		int vvanz=static_cast<int>(vanz);

		double xtol=pow(10,-5);

		int di=dinorm,i=dinorm;
		double xold=0;

		unordered_map<int,double> helpermap;
		Fmax=0;

		double dx=(vmax-vmin)/(vanz*dinorm);
		helpermap[0]=0.;
		Finv[0]=0.;
		for(double x=vmin+dx;x<=vmax;x=x+dx)
		{
			Fmax=Fmax+0.5*dx*(fpt(x-dx,rest...)+fpt(x,rest...));
			int ind=static_cast<int>(round((x-vmin)/dx));
			helpermap[ind]=Fmax;
		}
		unordered_map<int,double> *helpermap2=&helpermap;
		F=helpermap;
		double unc=pow(10,-4);
		double dF=unc*Fmax/vanz;

		while (i<dinorm*vanz)
		{
			int iold=i;
			double Ftest=static_cast<double>(i)*Fmax/(dinorm*vanz);
			double Ftest1=0;
			double Ftest2=Fmax;
			double Fakt=Ftest2;
			double xtest2=vmax;
			double xtest1=vmin+dx;

			while(abs(Ftest-Fakt)>dF and abs(xtest2-xtest1)>dx)
			{
				double xtmp2=0.5*(xtest1+xtest2);
				int xtmp2_approx=min(dinorm*vvanz,max(1,static_cast<int>(round((xtmp2-vmin)/dx))));
				double xx2;

				unordered_map<int,double> ::iterator it=helpermap2->find(xtmp2_approx);
				Fakt=it->second;

				if(Fakt>Ftest)
				{
					Ftest2=Fakt;
					xtest2=xtmp2;
				}
				else
				{
					Ftest1=Fakt;
					xtest1=xtmp2;
				}
			}
			Finv[i]=0.5*(xtest1+xtest2);

			if(Finv[i]-xold>=xtol)
			{
							xold=Finv[i];

				if(di>1){i=i-di; di=di/10;}
				
			}
			else{
							xold=Finv[i];
			if(i==(10*di)*floor(i/(10*di)) and di<dinorm){di=di*10;}

			}
			i=i+di;
		}
	}

	template<typename T, typename ...Rest> void fill_partitionfunction_v2(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,map<double,double> &Finv, double &Fmax,T x, Rest...rest )
	{
		double w=0;
		double dv=(vmax-vmin)/vanz;

		for(double v=vmin;v<=vmax;v=v+dv)
		{
			w=w+(dv/6.)*(f(v,rest...)+4.*f(v+0.5*dv,rest...)+f(v+dv,rest...));
			Finv[w]=v+dv;
		}
		Fmax=w;
	}