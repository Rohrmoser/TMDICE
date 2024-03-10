#include"deps.h"

#include "TMDICE_lib.h"
#include"globalconsts.h"

using namespace std;




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

unordered_map<int,double> discretize_grid(map<double,double>oldmap,double mapmax, int anz)
{
	unordered_map<int,double> newmap;
	map<double, double>::iterator it,it2;

	int anztmp=1,ii=0;
	newmap[0]=0.;
	double lastval,valmax;
	for (it = oldmap.begin(); it != oldmap.end(); it++)
	{
	it2=it;
	it2++;
	double ph1tmp=it->first,ph2tmp=it2->first;
	double z1tmp=it->second,z2tmp=it2->second;
	ii++;
	if(anz*ph1tmp/mapmax<anztmp and anz*ph2tmp/mapmax>anztmp )
	{
		lastval=z1tmp+(static_cast<double>(anztmp)*mapmax/static_cast<double>(anz)-ph1tmp)* (z2tmp-z1tmp)/(ph2tmp-ph1tmp);
		newmap[anztmp]=lastval;
		anztmp++;

	}
	if(ii>oldmap.size()-1){valmax=it->second;}
	}

	newmap[anztmp]=valmax;

	return newmap;
}

void fill_partitionfunction_v2( fkt f, double vmin,double vmax,double vanz,map<double,double> &Finv, double &Fmax)
{
    double w=0;
    double dv=(vmax-vmin)/vanz;
    for(double v=vmin;v<=vmax-dv;v=v+dv)
    {
		w=w+(dv/6.)*(f(v)+4.*f(v+0.5*dv)+f(v+dv));
 
    	Finv[w]=v+dv;
    }
    Fmax=w;
}

void fill_partitionfunction_v2( fkt f, double vmin,double vmax,double vanz,map<double,double> &Finv, map<double,double>&F,double &Fmax)
{
    double w=0;
    double dv=(vmax-vmin)/vanz;
    for(double v=vmin;v<=vmax-dv;v=v+dv)
    {
		w=w+(dv/6.)*(f(v)+4.*f(v+0.5*dv)+f(v+dv));
 
    	Finv[w]=v+dv;
		F[v+dv]=w;
    }
    Fmax=w;
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
		while(phinvtmp->find(static_cast<int>(rind))==phinvtmp->end() and rind>0 and fac<=dinorm /*and fac<=*//*anz*//*dinorm*/)
		{
			
			rind=floor(oldrind/fac)*fac;
			fac=fac*10;
		}
		if(phinvtmp->find(static_cast<int>(rind))!=phinvtmp->end()){xtmp1=phinvtmp->find(static_cast<int>(rind))->second;}
		double rind2=ceil(r);
		oldrind=rind2;
		fac=10.;
		while(phinvtmp->find(static_cast<int>(rind2))==phinvtmp->end() and rind2<=anzvle*dinorm and fac<=dinorm/*and fac<=dinorm*/)
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


double select(fkt2 Finv, double Fmax,double wmin,double wmax,double t)
{
	double r=Fmax*rand01();
	double w=Finv(r,t);
	w=max(wmin,w);
	w=min(wmax,w);
	return w;
}

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

void make_3c(map<double,double>f,double (&c1)[M],double (&c2)[M],double (&c3)[M],double a,double b1,double b2,double b)
{
	make_c(f,c1,a,b1,0.);
	make_c(f,c2,b1,b2,0.);
	make_c(f,c3,b2,b,0.);
}

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

void cheb_coeff_split(map<double,map<double,fkt2>>kernel,map<double,map<double,map<double,double>>>ph_max,map<double,map<double,map<double
,map<double,double>>>>phi_inv,int f1, int f2
, double ymin,double (&c1)[3][3][tanzmax][M],double (&c2)[3][3][tanzmax][M],double (&c3)[3][3][tanzmax][M]
,map<double,map<double,map<double,double>>> (&b1),map<double,map<double,map<double,double>>> (&b2),map<double,map<double,map<double,double>>> (&b),double t)
{
	double dt=(tmax-tmin)/tanz;
	{
		int tind=static_cast<int>((t-tmin)/dt);
			vector<double> bound;
			b[f1][f2][t]=ph_max[f1][f2].lower_bound(t)->second;
		if(b[f1][f2][t]>0)
		{
		double dk=pow(10,-6);
		b1[f1][f2][t]=b[f1][f2][t]*(0.5-dk);
		b2[f1][f2][t]=b[f1][f2][t]*(0.5+dk);
		double cc1[M],cc2[M],cc3[M];
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

double interpolate_cdf(double z, double z_max, unordered_map<int,double> cdf, int cdfanz,double cdfmax ,double z_min)
{
	double zind=((z-z_min)/(z_max-z_min))*cdfanz;
	double zmintmp=max(0.,floor(zind)),zmaxtmp=min(static_cast<double>(cdfanz),ceil(zind));
	double phmintmp=max(0.,cdf[static_cast<int>(zmintmp)]),phmaxtmp=min(cdfmax,cdf[static_cast<int>(zmaxtmp)]);
	return phmintmp+(zind-zmintmp)*(phmaxtmp-phmintmp);
}


void GetTransverseDirections(fourmom p1,fourmom &p2, fourmom &p3)
{
	p3.p0=0;
	p2.p0=0;

	fourmom e1(0,1,0,0),e2(0,0,1,0),e3(0,0,0,1);
	p1.p0=0;
	fourmom ptmp;
	ptmp=p1*(1./sqrt(p1.sk3(p1)));//unit vector;
	p1=ptmp;
	if(!(p1==e1) and !(p1==e2))
	{
		p2=e1+p1*(-e1.sk3(ptmp));
		p3=e2+p1*(-e2.sk3(ptmp));
	}
	else
	{
		if(p1==e1){p2=e2;p3=e3;}
		if(p1==e2){p2=e1;p3=e3;}
	}
	
}

fourmom Rotate(fourmom axis,fourmom vec)//vector is the vector that will be rotated;
{
	double pt=sqrt(pow(axis.px,2)+pow(axis.py,2)),p=sqrt(pt*pt+pow(axis.pz,2));
	double sinth=pt/p,costh=axis.pz/p;
	double sinph=axis.py/pt,cosph=axis.px/pt;
	if(pt==0){sinph=0.;cosph=1.;}
	fourmom vectemp(vec.p0,vec.px*costh+vec.pz*sinth,vec.py,-vec.px*sinth+vec.pz*costh);
	fourmom vectemp2(vectemp.p0,vectemp.px*cosph-vectemp.py*sinph,vectemp.px*sinph+vectemp.py*cosph,vectemp.pz);
	return vectemp2;
}