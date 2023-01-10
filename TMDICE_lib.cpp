#include"deps.h"

#include "TMDICE_lib.h"
#include"globalconsts.h"




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
			double nalt=n;
			n=300;
    double w=0;
    double dv=(vmax-vmin)/vanz;
    double fkt=0.001;

       for(double v=vmin;v<=/*=*/vmax-dv;v=v+dv)
    {
        double xx;
   //     w=w+integral(f,max(vmin,v-dv),v,xx);
		
			//w=w+0.5*dv*(f(v)+f(v+dv));
		w=w+(dv/6.)*(f(v)+4.*f(v+0.5*dv)+f(v+dv));
 
       Finv[w]=v+dv;
    }
    Fmax=w;
	n=nalt;
}


void fill_partitionfunction_v3( fkt f, double vmin,double vmax,double vanz,map<double,double> &Finv, double &Fmax)
{
			double nalt=n;
			n=1000;
    double w=0;
    double dv=(vmax-vmin)/vanz;
    double fkt=0.01;

	map<double,double> Finvhelp;

       for(double v=vmin;v<=/*=*/vmax;v=v+fkt*dv)
    {
        double xx;
 //       w=w+integral(f,max(vmin,v-dv),v,xx);
			w=w+0.5*dv*(f(v)+f(v+dv));
 
       Finvhelp[v]=w;
    }
    Fmax=w;
	double dw=Fmax/vanz;
	w=dw;
       for(double v=vmin;v<=/*=*/vmax;v=v+fkt*dv)
	{
		if(Finvhelp[v]<=w and Finvhelp[v]>w)
		{
			Finv[w]=v;
			w=w+dw;
		}
	}
	n=nalt;
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



/* vector<double> det_boundaries(fkt f, double y, double zmin,double zmax)
{
	vector<double> t;
	double dz=pow(10,-5);
	for(double z=zmin;z<zmax;z=z+dz )
	{
		if((f(z)>1./y and f(z+dz)<=1./y) or (f(z)<1./y and f(z+dz)>=1./y))
		{
			double zz;
			double ival=integral(f,zmin,z,zz);
			cout<<"ival="<<ival<<endl;
			t.push_back(ival);
		}
	}
	return t;
} */


void make_c(map<double,double> f,double (&c)[M],double a, double b,double eee)
{
	map<double,double>::iterator it2=f.end(),it_tmp;
	it2--;
	a=a+eee;
	b=b+eee;
	
	for(int j=0;j<=M-1;j++)
	{
		c[j]=0;
		for(double k=1.0;k<=m;k++)
		{
			double x=( 0.5*(a+b +cos(pi*(k-0.5)/m)*(b-a)) )-eee;
			it_tmp=f.lower_bound(x);
			x=it_tmp->first;
			c[j]=c[j]+((f[x]))*cos(pi*static_cast<double>(j)*(k-0.5)/m);
		}
		c[j]=c[j]*2.0/m;
	}

}

/* void make_3c(map<double,double>f,double (&c1)[M],double (&c2)[M],double (&c3)[M],double a,double b1,double b2,double b)
{
	make_c(f,c1,a,b1,0.);
	make_c(f,c2,b1,b2,0.);
	make_c(f,c3,b2,b,0.);
} */

double approx_f(double x, double a, double b,double (&c)[M],double eee)
{
	double y=((x+eee)-0.5*((a+eee)+(b+eee)))/(0.5*((b+eee)-(a+eee)));
	double f=0.5*c[0];
	for(int k=1;k<=M-1;k++)
	{
		f=f+c[k]*T(y,k);
	}
	return (f);
}

double interpolate_integer_inversegrid_1d(unordered_map<int,double> *phinvtmp,double xminus, double xplus, double r)
{
	double xtmp;
	double xtmp1=xminus,xtmp2=xplus;
	if(xminus==xplus)
	{
		return xminus;
	}
	else
	{
		double rind=floor(r), fac=10.,oldrind=rind;
		while(phinvtmp->find(static_cast<int>(rind))==phinvtmp->end() and rind>0 /*and fac<=*//*anz*//*dinorm*/)
		{
			
			rind=floor(oldrind/fac)*fac;
			fac=fac*10;
		}
		if(phinvtmp->find(static_cast<int>(rind))!=phinvtmp->end()){xtmp1=phinvtmp->find(static_cast<int>(rind))->second;}
		double rind2=ceil(r);
		oldrind=rind2;
		fac=10.;
		while(phinvtmp->find(static_cast<int>(rind2))==phinvtmp->end() and rind2<=anzvle*dinorm /*and fac<=dinorm*/)
		{
			
			rind2=ceil(oldrind/fac)*fac;
			fac=fac*10.;
		}
		if(phinvtmp->find(static_cast<int>(rind2))!=phinvtmp->end()){xtmp2=phinvtmp->find(static_cast<int>(rind2))->second;}
		xtmp2=min(xplus,xtmp2);
		xtmp1=max(xminus,xtmp1);
		if(rind2!=rind){xtmp=xtmp1+(r-static_cast<double>(rind))*(xtmp2-xtmp1)/(static_cast<double>(rind2-rind));}else{xtmp=xtmp1;}
		return xtmp;
	}
}

/*double select(fkt Finv, double Fmax,double wmin,double wmax)
{
	double r=Fmax*rand01();
	double w=Finv(r);
	w=max(wmin,w);
	w=min(wmax,w);
	return w;
}*/

double select(fkt2 Finv, double Fmax,double wmin,double wmax,double t)
{
	double r=Fmax*rand01();
	double w=Finv(r,t);
	w=max(wmin,w);
	w=min(wmax,w);
	return w;
}

//double c_t[3][3][tanzmax][M],c1_t[3][3][tanzmax][M],c2_t[3][3][tanzmax][M],c3_t[3][3][tanzmax][M];
//map<double,map<double,map<double,double>>>  b_t,b1_t,b2_t;


vector<double> det_boundaries(fkt f, double y, double zmin,double zmax)
{
	vector<double> t;
	double dz=pow(10,-6);
	for(double z=zmin;z<zmax;z=z+dz )
	{
		if((f(z)>1./y and f(z+dz)<=1./y) or (f(z)<1./y and f(z+dz)>=1./y))
		{
			double zz;
			double ival=integral(f,zmin,z,zz);
			t.push_back(ival);
		//	cout<<z<<" "<<ival<<endl;
		}
	}
	return t;
}

vector<double> det_boundaries(fkt2 f, double y, double zmin,double zmax,double t)
{
	vector<double> tt,tt2;
	double dz=pow(10,-5);
	for(double z=zmin;z<zmax;z=z+dz )
	{
//		cout<<z<<" "<<t<<" "<<f(z,t)<<" "<<f(z+dz,t)<<endl;
		if((f(z,t)>1./y and f(z+dz,t)<=1./y) or (f(z,t)<1./y and f(z+dz,t)>=1./y))
		{
			double zz;
			double ival=integral(f,zmin,z,zz,t);
			tt.push_back(ival);
			tt2.push_back(z);
		}
	}
	cout<<"			"<<tt[0]<<" "<<tt[1]<<endl;
	cout<<"			"<<tt2[0]<<" "<<f(tt2[0],t)<<endl;
	cout<<"			"<<tt2[1]<<" "<<f(tt2[1],t)<<endl;
	return tt;
}

/*void make_c(map<double,double> f,double (&c)[M],double a, double b,double eee)
{
	map<double,double>::iterator it2=f.end(),it_tmp;
	it2--;
	a=a+eee;
	b=b+eee;
	
	for(int j=0;j<=M-1;j++)
	{
		c[j]=0;
		for(double k=1.0;k<=m;k++)
		{
			double x=( 0.5*(a+b +cos(pi*(k-0.5)/m)*(b-a)) )-eee;
			it_tmp=f.lower_bound(x);
			x=it_tmp->first;
			c[j]=c[j]+((f[x]))*cos(pi*static_cast<double>(j)*(k-0.5)/m);
		}
		c[j]=c[j]*2.0/m;
	}
}*/

void make_3c(map<double,double>f,double (&c1)[M],double (&c2)[M],double (&c3)[M],double a,double b1,double b2,double b)
{
	make_c(f,c1,a,b1,0.);
	make_c(f,c2,b1,b2,0.);
	make_c(f,c3,b2,b,0.);
}

/*double approx_f(double x, double a, double b,double (&c)[M],double eee)
{
	double y=((x+eee)-0.5*((a+eee)+(b+eee)))/(0.5*((b+eee)-(a+eee)));
	double f=0.5*c[0];
	for(int k=1;k<=M-1;k++)
	{
		f=f+c[k]*T(y,k);
	}
	return (f);
}*/

double approx_f(double x,double *c1,double *c2,double *c3,double a,double b1,double b2,double b,int M)
{
	double aa,bb,f;
	if(a<b)
	{
		if(x<=b1)
		{
			aa=a,bb=b1;
			double y=((x)-0.5*((aa)+(bb)))/(0.5*((bb)-(aa)));

			f=0.5*(*c1);
			for(int k=1;k<=M-1;k++)
			{
				f=f+(*(c1+k))*T(y,k);
			}
		}
		if(x>b1 and x<=b2)
		{
			aa=b1;bb=b2;
			double y=((x)-0.5*((aa)+(bb)))/(0.5*((bb)-(aa)));
			f=0.5*(*c2);
			for(int k=1;k<=M-1;k++)
			{
				f=f+(*(c2+k))*T(y,k);
			}
		}
		if(x>b2)
		{
			aa=b2;bb=b;
			double y=((x)-0.5*((aa)+(bb)))/(0.5*((bb)-(aa)));
			f=0.5*(*c3);
			for(int k=1;k<=M-1;k++)
			{
				f=f+(*(c3+k))*T(y,k);
			}
		}
	}
	else
	{
		f=0;
	}

	return f;
}

void cheb_coeff_split(map<double,map<double,fkt2>>kernel,map<double,map<double,map<double,double>>>ph_max,map<double,map<double,map<double,map<double,double>>>>phi_inv,int f1, int f2
, double ymin,double (&c1)[3][3][tanzmax][M],double (&c2)[3][3][tanzmax][M],double (&c3)[3][3][tanzmax][M]
,map<double,map<double,map<double,double>>> (&b1),map<double,map<double,map<double,double>>> (&b2),map<double,map<double,map<double,double>>> (&b),double t)
{
	double dt=(tmax-tmin)/tanz;
//	for(double t=tmin;t<=tmax/sqrt(xmin);t=t+dt)
	{
		int tind=static_cast<int>((t-tmin)/dt);
			vector<double> bound;
			b[f1][f2][t]=ph_max[f1][f2].lower_bound(t)->second;
//		cout<<"tada "<<ph_max[f1][f2].lower_bound(t)->first<<endl;
		if(b[f1][f2][t]>0)
		{
	//	bound=det_boundaries(kernel[f1][f2],ymin/kernel[f1][f2](1.0-2.*xeps,t),1*xeps,1.-xeps,t);
		//b1gg=bound.at(0);
		//b2gg=bound.at(1);
		//bgg=phmax[2][2];
		double dk=pow(10,-6);
		b1[f1][f2][t]=b[f1][f2][t]*(0.5-dk);//bound.at(1);
		b2[f1][f2][t]=b[f1][f2][t]*(0.5+dk);//bound.at(0);
//		b1[f1][f2][t]=bound.at(0);
//		b2[f1][f2][t]=bound.at(1);
		double cc1[M],cc2[M],cc3[M];
//		cout<<(phi_inv[f1][f2].lower_bound(t)->first)<<" "<<(phi_inv[f1][f2].lower_bound(t)->second)<<endl;
	//	double aa1,aa2;
	//	aa1=phi_inv[f1][f2].lower_bound(t)->first;aa2=phi_inv[f1][f2].lower_bound(t)->second;
		make_3c(phi_inv[f1][f2].lower_bound(t)->second,cc1,cc2,cc3,0.,b1[f1][f2][t],b2[f1][f2][t],b[f1][f2][t]);	
		for(int i=0;i<M;i++){c1[f1][f2][tind][i]=cc1[i];c2[f1][f2][tind][i]=cc2[i];c3[f1][f2][tind][i]=cc3[i];}
		}
		else
		{
			for(int i=0;i<M;i++){c1[f1][f2][tind][i]=0;c2[f1][f2][tind][i]=0;c3[f1][f2][tind][i]=0;}
		}
		bound.clear();
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




