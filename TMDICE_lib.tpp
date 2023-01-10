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

	template<typename T, typename ...Rest> void fill_partitionfunction_v9(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,unordered_map<int,double> &F,unordered_map<int,double> &Finv, double &Fmax,T x, Rest...rest )
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

				if(di>1){i=i-di; di=di/10;/*dF=dF/10.;*/}
				
			}
			else{
							xold=Finv[i];
			if(i==(10*di)*floor(i/(10*di)) and di<dinorm){di=di*10;/*dF=dF*10.;*//*cout<<"Seas"<<endl;*/}

			}
			i=i+di;
		}
	}

	template<typename T, typename ...Rest> void fill_partitionfunction_v8(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,unordered_map<int,double> &Finv, double &Fmax,T x, Rest...rest )
	{
		//double nalt=n;
		//n=100;
//		double xx;
//		Fmax=integral(f,vmin,vmax,xx,rest...);

		double (*fpt)(T, Rest...)=(f);
		
		int vvanz=static_cast<int>(vanz);

//		double pmax=3.;
		double xtol=pow(10,-5);

		//int dinorm=pow(10,pmax),
		int di=dinorm,i=dinorm;
		double xold=0;

		unordered_map<int,double> helpermap;
		Fmax=0;



		double dx=(vmax-vmin)/(vanz*dinorm);
		helpermap[0]=0.;
		Finv[0]=0.;
		for(double x=vmin+dx;x<=vmax;x=x+dx)
		{
			//Fmax=Fmax+0.5*dx*(f(x-dx,rest...)+f(x,rest...));
			Fmax=Fmax+0.5*dx*(fpt(x-dx,rest...)+fpt(x,rest...));
			int ind=static_cast<int>(round((x-vmin)/dx));
			helpermap[ind]=Fmax;
//			cout<<x<<" "<<Fmax<<" "<<helpermap.find(ind)->first<<endl;
			/*if(x>=xnext)
			{
				int ind=static_cast<int>(vanztmp*Fmax);
				Finv[ind]=x;
				xnext=xnext+xtol;
				cout<<x<<" "<<ind<<" "<<xnext<<" "<<Finv[ind]<<" "<<Fmax<<endl;
			}*/
		}
		unordered_map<int,double> *helpermap2=&helpermap;
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
				//cout<<vmin<<" "<<xtmp2<<" ";
				int xtmp2_approx=min(dinorm*vvanz,max(1,static_cast<int>(round((xtmp2-vmin)/dx))));
//cout<<xtmp2<<" "<<xtmp2_approx<<" "<</*helpermap.find(xtmp2_approx)->second<<*/endl;
				double xx2;

				unordered_map<int,double> ::iterator it=helpermap2->find(xtmp2_approx);
				/*if(it==helpermap.end())
				{
					Fakt=integral(f,vmin,xtmp2,xx2,rest...);
					helpermap[xtmp2]=Fakt;
				}
				else
				{*/
					Fakt=it->second;
				//}
				
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
//				cout<<i<<" "<<xtest1<<" "<<xtest2<<" "<<Fakt<<" "<<Ftest<<endl;
				//cout<<xtest<<" "<<Ftest<<endl;
			}
			Finv[i]=0.5*(xtest1+xtest2);

			if(Finv[i]-xold>=xtol)
			{
							xold=Finv[i];

				if(di>1){i=i-di; di=di/10;/*dF=dF/10.;*/}
				
			}
			else{
							xold=Finv[i];
			if(i==(10*di)*floor(i/(10*di)) and di<dinorm){di=di*10;/*dF=dF*10.;*//*cout<<"Seas"<<endl;*/}

			}
			i=i+di;
		}
	//	n=nalt;
	}

template<typename T, typename ...Rest> void fill_partitionfunction_v7(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,unordered_map<int,double> &Finv, double &Fmax,T x, Rest...rest )
{
	double pmax=2.;
	double vanztmp=vanz*pow(10,pmax);
	double vvanz=pow(10,6),xtol=pow(10,-6);

	unordered_map<double,double> helpermap;
	double xnext=vmin+xtol;
	double dx=(vmax-vmin)/vvanz;
	Fmax=0;
	for(double x=vmin+dx;x<vmax;x=x+dx)
	{
		Fmax=Fmax+0.5*dx*(f(x-dx,rest...)+f(x,rest...));
		helpermap[x]=Fmax;
		/*if(x>=xnext)
		{
			int ind=static_cast<int>(vanztmp*Fmax);
			Finv[ind]=x;
			xnext=xnext+xtol;
//		cout<<x<<" "<<ind<<" "<<xnext<<" "<<Finv[ind]<<" "<<Fmax<<endl;
		}*/
	}
	for(double x=vmin+dx;x<vmax;x=x+dx)
	{
	//	Fmax=Fmax+0.5*dx*(f(x-dx,rest...)+f(x,rest...));
	//	helpermap[x]=Fmax;
		if(x>=xnext)
		{
			int ind=static_cast<int>(vanztmp*helpermap[x]/Fmax);
			Finv[ind]=x;
			xnext=xnext+xtol;
//		cout<<x<<" "<<ind<<" "<<xnext<<" "<<Finv[ind]<<" "<<Fmax<<endl;
		}
	}		
//	cout<<Fmax<<endl;

//	double xnext=vmin+xtol;
//	double dv=Fmax/vanz;
//	double vnext=dv;
//	double xc1,xc2,vc;
	/*for(double x=vmin+dx;x<=vmax;x=x+dx)
	{
		/*if(helpermap[x]>vnext and helpermap[x-dx]<=vnext)
		{
			double currentx=x-dx;// + (vnext-helpermap[x-dx])*dx/(helpermap[x]-helpermap[x-dx]);
			Finv[static_cast<int>(vnext*vanztmp/Fmax)]=currentx;
			vnext=vnext+dv;
			xnext=currentx+xtol;
		}


		if(helpermap[x]/(Fmax/vanztmp)>floor(helpermap[x-dx]/(Fmax/vanztmp)) and helpermap[x-dx]/(Fmax/vanztmp)<=floor(helpermap[x-dx]/(Fmax/vanztmp)) )
		{
			xc1=x-dx;
			xc2=x;
			vc=floor(helpermap[x-dx]/(Fmax/vanztmp))*(Fmax/vanztmp);
		}

		if(x-dx<=xnext and x>xnext)
		{
			Finv[static_cast<int>(vc*vanztmp/Fmax)]=xc2;//1+(vc-helpermap[xc1])*dx/(helpermap[xc2]-helpermap[xc1]);
			xnext=xc2+xtol;
		}*/
	//	if(/*helpermap[x]>vnext or */x>xnext)
	//	{
			/*int ind=static_cast<int>(vanztmp*helpermap[x]/Fmax);
			double xcur=x;
			Finv[ind]=xcur;*/
	//		int ind;
	//		double xcur;
			/*if(helpermap[x]>vnext)
			{
				ind=static_cast<int>(vanztmp*vnext/Fmax);
				xcur=x-dx+(vnext-helpermap[x-dx])*dx/(helpermap[x]-helpermap[x-dx]);
				vnext=vnext+dv;
			}*/
	//		if(x>xnext /*and helpermap[x]<=vnext*/)
	//		{
	//			double vcur=helpermap[x];//-dx]+(xnext-x+dx)*(helpermap[x]-helpermap[x-dx])/dx;
	//			ind=static_cast<int>(vanztmp*vcur/Fmax);
	//			xcur=xnext;
	//			xnext=x+xtol;

	//		}
	//		Finv[ind]=xcur;
		//}




	//}
}

	template<typename T, typename ...Rest> void fill_partitionfunction_v6(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,unordered_map<int,double> &Finv, double &Fmax,T x, Rest...rest )
	{
		//double nalt=n;
		//n=100;
		double xx;
		Fmax=integral(f,vmin,vmax,xx,rest...);

		double unc=0.00001;
		double dF=unc*Fmax/vanz;
		int vvanz=static_cast<int>(vanz);

		double pmax=3.;
		double xtol=pow(10,-3);

		int dinorm=pow(10,pmax),di=dinorm,i=dinorm;
		double xold=0;

		unordered_map<double,double> helpermap;
		while (i<dinorm*vanz)
		{
			int iold=i;
			double Ftest=static_cast<double>(i)*Fmax/(dinorm*vanz);
			double Ftest1=0;
			double Ftest2=Fmax;
			double Fakt=Ftest2;
			double xtest2=vmax;
			double xtest1=vmin;


			while(abs(Ftest-Fakt)>dF)
			{
				double xtmp2=0.5*(xtest1+xtest2);
				//cout<<vmin<<" "<<xtmp2<<" ";

				double xx2;

				unordered_map<double,double> ::iterator it=helpermap.find(xtmp2);
				if(it==helpermap.end())
				{
					Fakt=integral(f,vmin,xtmp2,xx2,rest...);
					helpermap[xtmp2]=Fakt;
				}
				else
				{
					Fakt=it->second;
				}
				
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
//				cout<<i<<" "<<xtest1<<" "<<xtest2<<" "<<Fakt<<" "<<Ftest<<endl;
				//cout<<xtest<<" "<<Ftest<<endl;
			}
			Finv[i]=0.5*(xtest1+xtest2);

			if(Finv[i]-xold>=xtol)
			{
							xold=Finv[i];

				if(di>1){i=i-di; di=di/10;/*dF=dF/10.;*/}
				
			}
			else{
							xold=Finv[i];
			if(i==(10*di)*floor(i/(10*di)) and di<dinorm){di=di*10;/*dF=dF*10.;*//*cout<<"Seas"<<endl;*/}

			}
			i=i+di;
		}
	//	n=nalt;
	}	

	template<typename T, typename ...Rest> void fill_partitionfunction_v5(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,unordered_map<int,double> &Finv, double &Fmax,T x, Rest...rest )
	{
		//double nalt=n;
		//n=100;
		double xx;
		Fmax=integral(f,vmin,vmax,xx,rest...);

		double unc=0.0001;
		double dF=unc*Fmax/vanz;
		int vvanz=static_cast<int>(vanz);
		for(int i=1;i<vvanz;i++)
		{
			double Ftest=static_cast<double>(i)*Fmax/vanz;
			double Ftest1=0;
			double Ftest2=Fmax;
			double Fakt=Ftest2;
			double xtest2=vmax;
			double xtest1=vmin;
			while(abs(Ftest-Fakt)>dF)
			{
				double xtmp2=0.5*(xtest1+xtest2);
			//	cout<<vmin<<" "<<xtmp2<<" ";

				double xx2;
				Fakt=integral(f,vmin,xtmp2,xx2,rest...);
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
				//cout<<xtest1<<" "<<xtest2<<" "<<Fakt<<" "<<Ftest<<endl;
				//cout<<xtest<<" "<<Ftest<<endl;
			}
			Finv[i]=0.5*(xtest1+xtest2);
		//	cout<<i<<" "<<Finv[i]<<" "<<Ftest<<" "<<Fmax<<"              Hi"<<endl;
		}
	//	n=nalt;
	}

	template<typename T, typename ...Rest> void fill_partitionfunction_v2(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,map<double,double> &Finv, double &Fmax,T x, Rest...rest )
	{
			double nalt=n;
			n=1000;
		double w=0;
		double dv=(vmax-vmin)/vanz;
		double fkt=0.001;

		   for(double v=vmin;v</*=*/vmax;v=v+dv)
		{
			double xx;
//			w=w+integral(f,max(vmin,v-dv),v,xx,rest...);
			w=w+0.5*dv*(f(v,rest...)+f(v+dv,rest...));
			Finv[w]=v;
		}
		Fmax=w;
		n=nalt;
	}

	template<typename T, typename ...Rest> void fill_partitionfunction_v3(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,unordered_map<double,double> &Finv, double &Fmax,T x, Rest...rest )
	{
			double nalt=n;
			n=1000;
    double w=0;
    double dv=(vmax-vmin)/vanz;
    double fkt=0.0001;

	map<double,double> Finvhelp;

       for(double v=vmin;v<=/*=*/vmax;v=v+fkt*dv)
    {
        double xx;
       // w=w+integral(f,max(vmin,v-dv),v,xx);
			w=w+0.5*dv*fkt*(f(v,rest...)+f(v+fkt*dv,rest...));
 //cout<<w<<" "<<v<<endl;
    //   Finvhelp[v]=w;
    }
    Fmax=w;
	double dw=Fmax/vanz;
	w=dw;
	double ww=0;
	for(double v=vmin;v<=/*=*/vmax;v=v+fkt*dv)
    {

			ww=ww+0.5*dv*fkt*(f(v,rest...)+f(v+fkt*dv,rest...));
 //	cout<<w<<" "<<v<<" "<<ww<<endl;
       if(ww>=w and ww<w+dw)
	   {
		   Finv[w]=v;
		   cout<<Fmax<<" "<<w<<" "<<v<<" "<<ww<<"                  Hi!"<<endl;
		   w=w+dw;
	   }
    }
	n=nalt;
}


	template<typename T, typename ...Rest> void fill_partitionfunction_v4(  double f(  T x,Rest...rest), double vmin,  double vmax, double vanz,unordered_map<int,double> &Finv, double &Fmax,T x, Rest...rest )
	{
			double nalt=n;
			n=1000;
    double w=0;
    double dv=(vmax-vmin)/vanz;
    double fkt=0.001;

	map<double,double> Finvhelp;

       for(double v=vmin;v<=/*=*/vmax;v=v+fkt*dv)
    {
        double xx;
       // w=w+integral(f,max(vmin,v-dv),v,xx);
			w=w+0.5*dv*fkt*(f(v,rest...)+f(v+fkt*dv,rest...));
 //cout<<w<<" "<<v<<endl;
    //   Finvhelp[v]=w;
    }
    Fmax=w;
	double dw=Fmax/vanz;
	w=dw;
	double ww=0;
	int ii=1;
	for(double v=vmin;v<=/*=*/vmax;v=v+fkt*dv)
    {

			ww=ww+0.5*dv*fkt*(f(v,rest...)+f(v+fkt*dv,rest...));
 //cout<<w<<" "<<v<<" "<<ww<<endl;
       if(ww>=w and ww<w+dw)
	   {
		   Finv[ii]=v;
//		   cout<<Fmax<<" "<<w<<" "<<v<<"                  Hi!"<<endl;
		   ii++;
		   w=ii*dw;
		   
	   }
    }
	n=nalt;
}
