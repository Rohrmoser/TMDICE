#include"deps.h"
#include"TMDICE_lib.h"
#include"decor.h"
#include "globalconsts.h"
#include "functions_lib.h"
#include "TMDICE_base.h"
#include"vle.h"

using namespace std;

bool nokinconstraints=false;

double ktmin;
double probggvle;
map<double,map<double,unordered_map<int,double>>> phiinvvle,phi_z_vle;

map<double,map<double,double>> phmaxvle;
map<double,double> phivle;

double pqq(double z){if(checkdla==1){return (alphas/(2.*pi))*cf/(1.-z);}else{return (alphas/(2.*pi))*cf*(1.+z*z)/(1.-z);}}
double pqg(double z){if(checkdla==1){return 0*1.0/z;}else{return 1*(alphas/(2.*pi))*(0.5)*(z*z+pow(1.-z,2));}}
double pgq(double z){return 0*1.0/z;}
double pgg(double z){if(checkdla==1){return 0.;}else{return 1*(alphas/(2.*pi))*ca*((1.-z)/z+z/(1.-z)+z*(1.0-z));}}
map<double,map<double,fkt>>P={{1.0, {{1.0,pqq},{2.0,pgq}}},{-1.0, {{-1.0,pqq},{2.0,pgq}}},{2.0, {{1.0,pqg},{2.0,pgg}}}};

int anzphiinvvle=pow(10,3);
int anzphivle=pow(10,3);

void obtaincdfsvle()
{
	for(double i=2.;i>=1.;i--)
	{
		for(double j=2;j>=1;j--)
		{
			map<double,map<double,string>>pnam={{1,{{1,"qq"},{2,"qg"}}},{2,{{1,"gq"},{2,"gg"}}}};
			cout<<"Calculating "<<pnam[i][j]<<" distributions:"<<endl<<endl;
			map<double,double>phtmp,phztmp;
			fill_partitionfunction_v2(P[i][j],1*xepsvle,1.0-1*xepsvle,anzvle,phtmp,phztmp,phmaxvle[i][j]);
			phiinvvle[i][j]=discretize_grid(phtmp,phmaxvle[i][j], anzphiinvvle);
			phi_z_vle[i][j]=discretize_grid(phztmp,1.0-1*xepsvle, anzphivle);
            phtmp.clear();
		}
	}
	phivle[1.]=phmaxvle[1.][1.];
	phivle[2.]=phmaxvle[2.][2.]+phmaxvle[2.][1.];
	probggvle=phmaxvle[2.][2.]/phivle[2.];
}

double phqqvle(double z){if(phmaxvle[1][1]>0){return interpolate_cdf(z,phmaxvle[1][1],phiinvvle[1][1],anzphiinvvle);}else{return 0.;}}
double phqgvle(double z){if(phmaxvle[2][1]>0){return interpolate_cdf(z,phmaxvle[2][1],phiinvvle[2][1],anzphiinvvle);}else{return 0.;}}
double phggvle(double z){if(phmaxvle[2][2]>0){return interpolate_cdf(z,phmaxvle[2][2],phiinvvle[2][2],anzphiinvvle);}else{return 0.;}}
double phgqvle(double z){return 0.;}

map<double,map<double,fkt>>phiinvvle1={{1.0, {{1.0,phqqvle},{2.0,phgqvle}}},{-1.0, {{-1.0,phqqvle},{2.0,phgqvle}}},{2.0, {{1.0,phqgvle},{2.0,phggvle}}}};

double phqqvlez(double z){if(phmaxvle[1][1]>0){return interpolate_cdf(z, 1.0-xepsvle,phi_z_vle[1][1],anzphivle,phmaxvle[1][1],xepsvle);}else{return 0.;}}
double phqgvlez(double z){if(phmaxvle[2][1]>0){return interpolate_cdf(z, 1.0-xepsvle,phi_z_vle[2][1],anzphivle,phmaxvle[2][1],xepsvle);}else{return 0.;}}
double phggvlez(double z){if(phmaxvle[2][2]>0){return interpolate_cdf(z, 1.0-xepsvle,phi_z_vle[2][2],anzphivle,phmaxvle[2][2],xepsvle);}else{return 0.;}}
double phgqvlez(double z){return 0.;}

map<double,map<double,fkt>>phi_z_vle1={{1.0, {{1.0,phqqvlez},{2.0,phgqvlez}}},{-1.0, {{-1.0,phqqvlez},{2.0,phgqvlez}}},{2.0, {{1.0,phqgvlez},{2.0,phggvlez}}}};

bool stop_with_ktmin(TMDICEparticle p)
{
	bool cond;
	cond=(p.kt12()<=ktmin and p.adr!="1");;
	
	return cond;
}

bool stop_when_tf_exceeds_tmax(TMDICEparticle p)
{
	return (p.tf()>=tmax);
}

bool stop_when_tf_exceeds_td(TMDICEparticle p)
{
	bool cond;
	cond=(p.tf()>=p.td());
	return cond; 
}

bool stop_when_td_exceeds_tmax(TMDICEparticle p)
{
	return (p.td()>=tmax);
}

bool stop_non_resolved_but_resolvable(TMDICEparticle p)
{
	bool cond1,cond2,cond3;
	cond1=stop_with_ktmin(p);
	cond2=stop_when_tf_exceeds_td(p);
	cond3=stop_when_td_exceeds_tmax(p);
	return (cond1 or cond2 or cond3);
}

bool stop_with_ktmin(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	bool cond;
	cond=(p1.kt12()<=ktmin or p2.kt12()<=ktmin );
	
	return cond;
}

bool stop_when_tf_exceeds_tmax(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	return(p1.tf()>=tmax);
}

bool stop_when_tf_exceeds_td(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	bool cond;
	cond=(p1.tf()>=p1.td());
	return cond;
}

bool stop_when_td_exceeds_tmax(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	return (p1.td()>=tmax) ;
}

bool stop_non_resolved_but_resolvable(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	bool cond1,cond2,cond3;
	cond1=stop_with_ktmin(p,p1,p2);
	cond2=stop_when_tf_exceeds_td(p,p1,p2);
	cond3=stop_when_td_exceeds_tmax(p,p1,p2);
	return (cond1 or cond2 or cond3 or p1.Q==0);
}

bool stop_with_qmin(TMDICEparticle p)
{
	return (p.Q<Q2min);
}

bool ltxminvle(TMDICEparticle p)
{
	return (p.x<=xminvle);
}

void settoq(TMDICEparticle &p,double q)
{
	p.Q=q;
}

//VLEevent:
VLEevent::VLEevent()
{
	setevolmin(Q2min);
	setsetpval(settoq);
	setstopcond(stop_with_qmin);
	setdumpcond(ltxminvle);
}

VLEevent::VLEevent( fktpartbool stopcondtmp, fktpartbool3 nosplitcondtmp,bool filltempstmp,fktpartvoid pvaltmp, fktpartbool dumpcondtmp,double minq)
{
	setevolmin(minq);
	setsetpval(pvaltmp);
	setstopcond(stopcondtmp);
	setdumpcond(dumpcondtmp);
	filltemps=filltempstmp;	
	setnosplitcond(nosplitcondtmp);
}

double xaobound(double qtmp,double thold,double x)
{
	double ee=x*emax;
	double xboundtmp;
	if(1.-2.0*qtmp*qtmp/(ee*ee*(1.0-cos(thold)))>0.){xboundtmp=0.5*sqrt(1.-2.0*qtmp*qtmp/(ee*ee*(1.0-cos(thold))));}else{xboundtmp=0.;}

	return min(0.5-xepsvle,xboundtmp);
}

void vleboundaries(double q, double thold,double x, int & ixminus, int& ixplus, double & xminus, double & xplus)
{
	double dxbound;
	if(ao_mode==0){dxbound=xaobound(q,R,x);}
	if(ao_mode==1){dxbound=xaobound(q,thold,x);}

	xplus=0.5+dxbound;
	xminus=0.5-dxbound;

	double dx=(1.0-2*xepsvle)/(anzvle*static_cast<double>(dinorm));
	
	ixplus=min(static_cast<int>(anzvle*static_cast<double>(dinorm))-1,static_cast<int>(round((xplus-xepsvle)/dx)));
	ixminus=max(0,static_cast<int>(round((xminus-xepsvle)/dx)));
}

double VLEevent::evol_select(TMDICEparticle p)
{
	//cout<<"qsel"<<flush<<endl;
	TMDICEparticle* pold=p.old_part;
	double f=p.typ;
	double r1,r2,wwww;
	double z=p.x/(pold->x);
	double t1=pow(min(p.x*emax,pold->Q),2);
	double t2;
	double phi1=phivle[abs(f)];
	if(ao_mode!=-1)
	{
		double theta=p.th_old;
			do
			{
				
				r1=rand01();
				t2=exp(log(t1)+(1./phi1)*log(r1));

				double ttmp=t2;

				int ixminus, ixplus;
				double xminus, xplus;
				
				vleboundaries(  sqrt(ttmp),  theta, p.x,  ixminus, ixplus,xminus, xplus);
				

				double phim2=phi_z_vle1[abs(f)][1.](xplus)+phi_z_vle1[abs(f)][2.](xplus);
				double phim1=phi_z_vle1[abs(f)][1.](xminus)+phi_z_vle1[abs(f)][2.](xminus);

				double phim=phim2-phim1;
				if(xplus==xminus or ixplus==ixminus){phim=0; }		
				wwww=phim/phi1;
				t1=t2;
				r2=rand01();
			}while(wwww<r2 and t2>Q2min*Q2min);
	}
	if(ao_mode==-1)
	{
		double r=rand01(),tfin;
		t2=exp(log(t1)+(1./phi1)*log(r));
	}
	return sqrt(t2);
}

double VLEevent::typselect(TMDICEparticle p)
{
	double f1=p.typ;

	double f2;
	double r=rand01();
	if(abs(f1)==1.0){f2=f1;}

	if(f1==2.0)
	{
		double rcp;
		rcp=probggvle;
		if(ao_mode!=-1)
		{
			int ixminus, ixplus;
			double xminus, xplus;		
			vleboundaries( p.Q,  p.th_old, p.x,  ixminus, ixplus,xminus, xplus);

			double phim2=phi_z_vle1[2.][1.](xplus)+phi_z_vle1[2.][2.](xplus);
			double phim1=phi_z_vle1[2.][1.](xminus)+phi_z_vle1[2.][2.](xminus);

			double phimgg2=phi_z_vle1[2.][2.](xplus);
			double phimgg1=phi_z_vle1[2.][2.](xminus);		
			rcp=(phimgg2-phimgg1)/(phim2-phim1);	
			if(xplus==xminus or ixplus==ixminus ){rcp=1;}		

		}

		if(r<rcp){f2=2.0;}else{double r2=rand01(); if(r2<0.5){f2=1.0;}else{f2=-1.0;} }
	}
	return f2;
}

double VLEevent::xselect(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
//	cout<<"xsel"<<flush<<endl;
	double f1=p.typ;
	double f2=p1.typ;

	double xtmp;
	double r0,rmax;

	double xminus=xepsvle, xplus=1.0-xepsvle;
	if(ao_mode==-1)
	{
		r0=0;
		rmax=anzvle*static_cast<double>(dinorm);
		rmax=phmaxvle[abs(f1)][abs(f2)];
	}
	if(ao_mode!=-1)
	{
		int ixminus, ixplus,ixintminus,ixintplus;
		vleboundaries( p.Q,  p.th_old, p.x,  ixminus, ixplus,xminus, xplus);

		double phim2=phi_z_vle1[abs(f1)][abs(f2)](xplus);
		double phim1=phi_z_vle1[abs(f1)][abs(f2)](xminus);

		
		double phi1=phmaxvle[abs(f1)][abs(f2)];
	
		double rtmp=1;

		r0=rtmp*(phim1);
		rmax=rtmp*((phim2-phim1));

	}	
	double r=r0+rand01()*rmax;

	if(xplus==xminus){xtmp=0.5;}		
	else{
		xtmp=phiinvvle1[abs(f1)][abs(f2)](r);
		}
	return xtmp;
}

void RotateToP(TMDICEparticle p, TMDICEparticle &p1, double phirel)
{
	p1.p.p0=p1.x*emax;
	p1.p.px=p1.ktrel*cos(phirel);
	p1.p.py=p1.ktrel*sin(phirel);
	p1.p.pz=sqrt(p1.p.p0*p1.p.p0-p1.Q*p1.Q-p1.ktrel*p1.ktrel);


	double thtmp=acos(p.p.pz/sqrt(p.p.p0*p.p.p0-p.Q*p.Q)),phtmp=atan2(p.p.py,p.p.px);
	double p1ztmp=p1.p.pz*cos(thtmp)-p1.p.px*sin(thtmp);
	double p1xtmp=p1.p.pz*sin(thtmp)+p1.p.px*cos(thtmp);
	double p1ytmp=p1.p.py;
	p1.p.px=p1xtmp*cos(phtmp)-p1ytmp*sin(phtmp);
	p1.p.py=p1xtmp*sin(phtmp)+p1ytmp*cos(phtmp);
	p1.p.pz=p1ztmp;
}

void VLEevent::vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2)
{
	//cout<<"vertex"<<flush<<endl;
	p1=p;
	p2=p;

	p1.old_part=&p;
	p2.old_part=&p;
	

	p1.typ=typselect(p);	
	p2.typ=typ2[p.typ][p1.typ];
	double kttmp2,z;
	do{
		p1.x=xselect(p,p1,p2);	
		if(checkdla==0)
		{
			p2.x=1.0-p1.x;//this is now only for selection of z, i.e.: relative fraction
			p1.x=p1.x*p.x;//this is for total energy fraction x
			p2.x=p2.x*p.x;
		}

		if(checkdla==1)
		{p2.x=1.-p1.x;p1.x=1.;}
		determineadress(p1,p2);

		p1.Q=evol_select(p1);
		p2.Q=evol_select(p2);

		z=p1.x/p.x;
		kttmp2=p.Q*p.Q*z*(1.-z)-(1.-z)*p1.Q*p1.Q-z*p2.Q*p2.Q;
		kttmp2=kttmp2-0.25*(pow(p.Q,4.)+pow(p1.Q,4)+pow(p2.Q,4)-2.*pow(p1.Q*p.Q,2)-2.*pow(p.Q*p2.Q,2)-2.*pow(p1.Q*p2.Q,2))/pow(p.x*emax,2);
		kttmp2=kttmp2/(1.-pow(p.Q/(p.x*emax),2.));
	}while(kttmp2<0 and checkdla!=1 and nokinconstraints==false);

	p1.ktrel=sqrt(kttmp2);
	p2.ktrel=sqrt(kttmp2);

	p1.kt=p.kt*p1.x;
	p2.kt=p.kt*p2.x;

	p1.th_old=p2.th12();
	p2.th_old=p2.th12();

	double phirel=rand01()*2.*pi;	

	if(rotate==1)
	{
		p1.p.p0=p1.x*emax;
		p1.p.px=p1.ktrel*cos(phirel);
		p1.p.py=p1.ktrel*sin(phirel);
		p1.p.pz=sqrt(p1.p.p0*p1.p.p0-p1.Q*p1.Q-p1.p.px*p1.p.px-p1.p.py*p1.p.py);

		p2.p.p0=p2.x*emax;
		p2.p.px=p2.ktrel*cos(phirel+pi);
		p2.p.py=p2.ktrel*sin(phirel+pi);
		p2.p.pz=sqrt(p2.p.p0*p2.p.p0-p2.Q*p2.Q-p2.p.px*p2.p.px-p2.p.py*p2.p.py);

		fourmom tmp1=Rotate(p.p,p1.p);
		p1.p=tmp1;

		fourmom tmp2=Rotate(p.p,p2.p);
		p2.p=tmp2;
	}
	else
	{
		p1.p.p0=p1.x*emax;
		p1.p.px=p1.ktrel*cos(phirel);
		p1.p.py=p1.ktrel*sin(phirel);
		p1.p.pz=sqrt(p1.p.p0*p1.p.p0-p1.Q*p1.Q-p1.p.px*p1.p.px-p1.p.py*p1.p.py);

		p2.p.p0=p2.x*emax;
		p2.p.px=p2.ktrel*cos(phirel+pi);
		p2.p.py=p2.ktrel*sin(phirel+pi);
		p2.p.pz=sqrt(p2.p.p0*p2.p.p0-p2.Q*p2.Q-p2.p.px*p2.p.px-p2.p.py*p2.p.py);
		
	}
	p1.kt=p1.pt();
	p2.kt=p2.pt();

	p1.phik=p1.phi();
	p2.phik=p1.phi();	
}