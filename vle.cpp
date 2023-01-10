
#include"deps.h"
#include"TMDICE_lib.h"
#include"decor.h"
#include "globalconsts.h"
#include "functions_lib.h"
#include "TMDICE_base.h"
#include"vle.h"

double cqq_dglap[M], cqg_dglap[M],cgg_dglap[M];
double cqq1_dglap[M], cqg1_dglap[M],cgg1_dglap[M];
double cqq2_dglap[M], cqg2_dglap[M],cgg2_dglap[M];
double cqq3_dglap[M], cqg3_dglap[M],cgg3_dglap[M];
double b1gg_dglap,b2gg_dglap,bgg_dglap;
double b1qg_dglap,b2qg_dglap,bqg_dglap;
double b1qq_dglap,b2qq_dglap,bqq_dglap;

double ktmin;
double probggvle;
map<double,map<double,unordered_map<int,double>>> phiinvvle,phi_z_vle;

map<double,map<double,double>> phmaxvle;
map<double,double> phivle;

double pqq(double z){if(checkdla==1){return (alphas/(2.*pi))*cf/(1.-z);}else{return (alphas/(2.*pi))*cf*(1.+z*z)/(1.-z);}}
double pqg(double z){if(checkdla==1){return 0*1.0/z;}else{return (alphas/(2.*pi))*(0.5)*(z*z+pow(1.-z,2));}}
double pgq(double z){return 0*1.0/z;}
double pgg(double z){if(checkdla==1){return 0.;}else{return (alphas/(2.*pi))*ca*((1.-z)/z+z/(1.-z)+z*(1.0-z));}}
map<double,map<double,fkt>>P={{1.0, {{1.0,pqq},{2.0,pgq}}},{-1.0, {{-1.0,pqq},{2.0,pgq}}},{2.0, {{1.0,pqg},{2.0,pgg}}}};



void setVLE()
{
	disclaimer();

	cout<<endl<<"Calculating partition functions: This may take a while."<<endl;
	gevfm=1.0/0.19732;	
	xepsvle=1.*pow(10,-4);
	xminvle=1.*pow(10,-4);
	anzvle=1.*pow(10,5);
	cout<<flush<<anzvle<<" "<<dinorm<<flush<<endl;
	kt0=0.;
	x0=1.;
	typ0=2;
	N_F=3.;
	ao_mode=1;
	R=pi;
	checkdla=0;
	setbias=0;
	simplifythkt=0;

	map<string,double*> invars2;

	invars2={{"xepsvle",&xepsvle},{"xminvle",&xminvle},{"kt1",&kt0},{"x1",&x0},{"typ1",&typ0},{"nc",&nc},
	{"alphabar",&alphabar},{"alphas",&alphas},{"emax",&emax},{"theta0orR",&theta0orR},{"R",&R},{"Q2max",&Q2max},{"Q2min",&Q2min},{"qhat",&qhat},{"tmin",&tmin},{"tmax",&tmax}
	,{"ao_mode",&ao_mode},{"ktmin",&ktmin},{"checkdla",&checkdla},{"simplifythkt",&simplifythkt},{"setbias",&setbias},{"p1x",&p1x},{"p1y",&p1y},{"p1z",&p1z}};

	setparams(invars2);

	ca=nc;
	cf=(nc*nc-1.)/(2.*nc);

	colorfact=pow(cf/ca,1.0);

	if(invars.find("alphabar")==invars.end())
	{
		alphabar=invars["alphas"]/pi;
	}
	else
	{
		alphas=alphabar*pi;
	}
	qhatgev=qhat/gevfm;//inGeV^3
	cout<<"alphas="<<alphas<<endl;

	map<double,map<double,double>> ymap={{1,{{1,500000}}},{2,{{1,1000},{2,1*pow(10,10)}}}};
	//			map<double,map<double,map<double,double>>> phiinv,difffkt;

	for(double i=2.;i>=1.;i--)
	{
		for(double j=2;j>=1;j--)
		{
			map<double,map<double,string>>pnam={{1,{{1,"qq"},{2,"qg"}}},{2,{{1,"gq"},{2,"gg"}}}};
			cout<<"Calculating "<<pnam[i][j]<<" distributions:"<<endl<<endl;
			double xxx;
			fill_partitionfunction_v9(P[i][j],1*xepsvle,1.0-1*xepsvle,anzvle,phi_z_vle[i][j],phiinvvle[i][j],phmaxvle[i][j],xxx);
			//fill_partitionfunction_v2(P[i][j],1*xepsvle,1.0-1*xepsvle,anzvle,phiinvvle[i][j],phmaxvle[i][j]);

		}
	}
	phivle[1.]=phmaxvle[1.][1.];
	phivle[2.]=phmaxvle[2.][2.]+phmaxvle[2.][1.];
	probggvle=phmaxvle[2.][2.]/phivle[2.];

	srand (time(NULL));
	cout<<"Initialization finished"<<endl;
}

bool stop_with_ktmin(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	bool cond;
	if(setbias==1)
	{
		cond=(p.Q<=2.*ktmin);
	}
	else
	{
		cond=(p1.kt12()<=ktmin or p2.kt12()<=ktmin /*and p.adr!="1"*/);
	}
	return cond;
}

bool stop_with_ktmin(TMDICEparticle p)
{
	bool cond;
	if(setbias==1)
	{
		cond=(p.Q<=2.*ktmin);
	}
	else
	{
		cond=(p.kt12()<=ktmin and p.adr!="1");;
	}
	return cond;
}

bool stop_when_tf_exceeds_tmax(TMDICEparticle p)
{
	return (p.tf()>=tmax);//and(p.adr!="1");;
}

bool stop_when_tf_exceeds_tmax_and_ktmin(TMDICEparticle p)
{
	bool cond;
	if (setbias==0){cond=(p.tf()>=tmax)or ((p.adr!="1") and p.kt12()<=ktmin);}
	if (setbias==1){cond=(p.tf()>=tmax)or ((p.adr!="1") and p.Q<=2.*ktmin);}
	return cond;
}

bool stop_when_td_exceeds_tmax_and_ktmin(TMDICEparticle p)
{
	bool cond;
	if (setbias==0){cond=(p.td()>=tmax)or ((p.adr!="1") and p.kt12()<=ktmin);}
		double thetac=pow(12./(qhatgev*pow(tmax*gevfm,3)),0.5);
	if (setbias==1){cond=(p.td()>=tmax)or ((p.adr!="1") and p.kt12()<=2.*ktmin)or p.Q<=pow(ktmin*p.x*emax*thetac,0.5) ;}
	return cond;
}

bool stop_when_tf_exceeds_td(TMDICEparticle p)
{
	bool cond;
	if(setbias==0){cond=(p.tf()>=p.td());}//and(p.adr!="1");;
	if(setbias==1){cond=(p.Q<=pow(8.*p.x*emax*qhatgev/3.,0.25));}
	return cond; 
}

bool stop_when_td_exceeds_tmax(TMDICEparticle p)
{
	return (p.td()>=tmax);//and(p.adr!="1");;
}

bool stop_non_resolved_but_resolvable(TMDICEparticle p)
{
	bool cond1,cond2,cond3;
	cond1=stop_with_ktmin(p);
	cond2=stop_when_tf_exceeds_td(p);
	cond3=stop_when_td_exceeds_tmax(p);
	return (cond1 or cond2 or cond3);
}

bool stop_when_tf_exceeds_tmax(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	return(p1.tf()>=tmax);//and(p.adr!="1");;
}

bool stop_when_tf_exceeds_tmax_and_ktmin(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
//	return (p1.tf()>=tmax)or (p2.tf()>=tmax)or (p1.kt12()<=ktmin/* and(p.adr!="1")*/);;
	bool cond;
	if (setbias==0){cond=(p1.tf()>=tmax)or ( p1.kt12()<=ktmin);}
	if (setbias==1){cond=(p1.tf()>=tmax)or ( p.Q<=2.*ktmin);}
	return cond;
}

bool stop_when_tf_exceeds_td(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
//	return (p1.tf()>=p1.td());//and(p.adr!="1");;
	bool cond;
	if(setbias==0){cond=(p1.tf()>=p1.td());}//and(p.adr!="1");;
	if(setbias==1){cond=(p.Q<=pow(8.*p.x*emax*qhatgev/3.,0.25));}
	return cond;
}

bool stop_when_td_exceeds_tmax(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	return (p1.td()>=tmax) ;//and(p.adr!="1");;
}

bool stop_when_td_exceeds_tmax_and_ktmin(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	bool cond;
	if (setbias==0){cond=(p1.td()>=tmax)or ( p1.kt12()<=ktmin);}
	double z=p1.x/p.x;
	double thetac=pow(12./(qhatgev*pow(tmax*gevfm,3)),0.5);
	if (setbias==1){cond=(p1.td()>=tmax)or ( p.Q*sqrt(z*(1.-z))<=ktmin)or p.Q<=pow(ktmin*p.x*emax*thetac,0.5) ;}
	return cond;
}

bool stop_non_resolved_but_resolvable(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	bool cond1,cond2,cond3;
	cond1=stop_with_ktmin(p,p1,p2);
	cond2=stop_when_tf_exceeds_td(p,p1,p2);
	cond3=stop_when_td_exceeds_tmax(p,p1,p2);
	return (cond1 or cond2 or cond3);
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

double bias_ktmin(TMDICEparticle p)
{
	return 0.5*sqrt(max(0.,1.-pow(2.*ktmin/p.Q,2)));
}

double bias_tftd(TMDICEparticle p)
{
	return 0.5*sqrt(max(0.,1.-8.*p.x*emax*qhatgev/(2.*pow(p.Q,4.))));
}

double bias_tftd_ktmin(TMDICEparticle p)
{
	return min(bias_tftd(p),bias_ktmin(p));
}

double bias_tdtl(TMDICEparticle p)
{
	double val= -0.5*sqrt(max(0.,1.0-p.Q*p.Q*qhatgev*pow(tmax*gevfm,3)/(3*pow(p.x*emax,2))) );
	cout<<val<<endl;
	return val;
}

//VLEevent:
VLEevent::VLEevent()
{
	setevolmin(Q2min);
	setsetpval(settoq);
	setstopcond(stop_with_qmin);
	setdumpcond(ltxminvle);
	setselectionbias(no_bias);
}

VLEevent::VLEevent( fktpartbool stopcondtmp, fktpartbool3 nosplitcondtmp,fktpart selectionbiastmp,fktpartvoid pvaltmp, fktpartbool dumpcondtmp,double minq)
{
	setevolmin(minq);
	setsetpval(pvaltmp);
	setstopcond(stopcondtmp);
	setdumpcond(dumpcondtmp);	
	setnosplitcond(nosplitcondtmp);
	setselectionbias(selectionbiastmp);
}

double xaobound(double qtmp,double thold,double x)
{
	double ee=x*emax;
	double xboundtmp;
	if(1.-2.0*qtmp*qtmp/(ee*ee*(1.0-cos(thold)))>0.){xboundtmp=0.5*sqrt(1.-2.0*qtmp*qtmp/(ee*ee*(1.0-cos(thold))));}else{xboundtmp=0.;}//0.5-xepsvle;}
	return min(0.5-xepsvle,xboundtmp);
}

void VLEevent::vleboundaries(double q, double thold,double x, int & ixminus, int& ixplus, double & xminus, double & xplus,
 int & ixintminus, int& ixintplus, double & xintminus, double & xintplus)
{
	double dxbound;
	if(ao_mode==0){dxbound=xaobound(q,R,x);}
	if(ao_mode==1){dxbound=xaobound(q,thold,x);}

	TMDICEparticle p1;
	p1.Q=q;
	p1.x=x;
	double detx=selectionbias(p1);
	xplus=0.5+dxbound;
	xminus=0.5-dxbound;
	xintminus=0.5;
	xintplus=0.5;
//	cout<<xminus<<" "<<xintminus<<" "<<xintplus<<" "<<xplus<<endl;

	if(detx>=0)
	{
		xplus=min(xplus,1.-detx);
		xminus=max(xminus,detx);
	}
	else
	{
		xintminus=max(xminus,min(xintminus,0.5+detx));
		xintplus=min(xplus,max(0.5,0.5-detx));
	}
//	cout<<xminus<<" "<<xintminus<<" "<<xintplus<<" "<<xplus<<endl;
//	cout<<thold<<endl;
	//cout<<"ao_mode="<<ao_mode<<endl;

	double dx=(1.0-2*xepsvle)/(anzvle*static_cast<double>(dinorm));
	
	ixplus=min(static_cast<int>(anzvle*static_cast<double>(dinorm))-1,static_cast<int>(round((xplus-xepsvle)/dx)));
	ixminus=max(0,static_cast<int>(round((xminus-xepsvle)/dx)));

	ixintminus=max(0,static_cast<int>(round((xintminus-xepsvle)/dx)));
	ixintplus=min(static_cast<int>(anzvle*static_cast<double>(dinorm))-1,static_cast<int>(round((xintplus-xepsvle)/dx)));
}

void vleboundaries1(double q, double thold,double x, int & ixminus, int& ixplus, double & xminus, double & xplus)
{
	double dxbound;
	if(ao_mode==0){dxbound=xaobound(q,R,x);}
	if(ao_mode==1){dxbound=xaobound(q,thold,x);}

	//TMDICEparticle p1;
	//p1.Q=q;
	//dxbound=min(dxbound,selectionbias(p1));
	xplus=0.5+dxbound;
	xminus=0.5-dxbound;

	double dx=(1.0-2*xepsvle)/(anzvle*static_cast<double>(dinorm));
	
	ixplus=min(static_cast<int>(anzvle*static_cast<double>(dinorm))-1,static_cast<int>(round((xplus-xepsvle)/dx)));
	ixminus=max(0,static_cast<int>(round((xminus-xepsvle)/dx)));
}

double VLEevent::evol_select(TMDICEparticle p)
{
//	cout<<"qsel"<<flush<<endl;
	double f=p.typ;
	double colorfact2=1.;
	if(abs(f)==1){colorfact2=colorfact;}
	double r1,r2,wwww;
	double z=p.x/p.x_old;
	double t1=min(pow(p.x*emax,2),p.Q_old*p.Q_old);
	//double t1=pow(p.Q_old,2);
	double t2;
	double phi1=phivle[abs(f)];
	if(ao_mode!=-1)
	{
		do
		{
			
			r1=rand01();
			t2=exp(log(t1)+(1./phi1)*log(r1));

			double ttmp=t2;
//			cout<<"ttmp="<<(ttmp)<<flush<<endl;

			int ixminus, ixplus,ixintminus,ixintplus;
			double xminus, xplus,xintminus,xintplus;
			
			vleboundaries(  sqrt(ttmp),  p.th_old, p.x,  ixminus, ixplus,xminus, xplus,ixintminus,ixintplus,xintminus,xintplus);
			
			unordered_map<int,double> *phimtmp1=&phi_z_vle[abs(f)][1.],*phimtmp2=&phi_z_vle[abs(f)][2.];
	//		cout<<xminus<<" "<<xplus<<"		"<<ixminus<<" "<<ixplus<<" "<<anzvle*dinorm<<" "<<phimtmp2->size()<<flush<<endl;
			double phim2=(*phimtmp1).find(ixplus)->second+(*phimtmp2).find(ixplus)->second;
			double phim1=(*phimtmp1).find(ixminus)->second+(*phimtmp2).find(ixminus)->second;

			double phim4=(*phimtmp1).find(ixintplus)->second+(*phimtmp2).find(ixintplus)->second;
			double phim3=(*phimtmp1).find(ixintminus)->second+(*phimtmp2).find(ixintminus)->second;

			double phim=phim2-/*(phim4-phim3)-*/phim1;
			if(xplus==xminus or ixplus==ixminus/* or (ixintminus==ixminus and ixintplus==ixplus)*/){phim=0; /*cout<<"hi"<<flush<<endl;*/}		
			wwww=phim/phi1;
			t1=t2;
			r2=rand01();
		}while(wwww<r2 and t2>Q2min*Q2min);//check later!!!
	}
	if(ao_mode==-1)
	{
		double r=rand01(),tfin;
		t2=exp(log(t1)+(1./phi1)*log(r));
	}
//	cout<<"t="<<(t2)<<flush<<endl;
	return sqrt(t2);
}

double VLEevent::typselect(TMDICEparticle p)
{
//	cout<<"typsel"<<flush<<endl;
	double f1=p.typ;

	double f2;
	double r=rand01();
	if(f1==1.0){f2=1.0;}
	if(f1==-1.0){f2=-1.0;}

	if(f1==2.0)
	{
		double rcp;
		rcp=probggvle;
		if(ao_mode!=-1)
		{
			unordered_map<int,double> *phimtmp1=&phi_z_vle[2.][1.],*phimtmp2=&phi_z_vle[2.][2.];	

			int ixminus, ixplus,ixintminus,ixintplus;
			double xminus, xplus,xintminus,xintplus;		
			vleboundaries( p.Q,  p.th_old, p.x,  ixminus, ixplus,xminus, xplus,ixintminus,ixintplus,xintminus,xintplus);

			rcp=((*phimtmp2).find(ixplus)->second-(*phimtmp2).find(ixintplus)->second+(*phimtmp2).find(ixintminus)->second-(*phimtmp2).find(ixminus)->second)/
			((*phimtmp2).find(ixplus)->second+(*phimtmp1).find(ixplus)->second
			//-(*phimtmp2).find(ixintplus)->second-(*phimtmp1).find(ixintplus)->second
			//+(*phimtmp2).find(ixintminus)->second+(*phimtmp1).find(ixintminus)->second
			-(*phimtmp2).find(ixminus)->second-(*phimtmp1).find(ixminus)->second);
			if(xplus==xminus or ixplus==ixminus /*or (ixintminus==ixminus and ixintplus==ixplus)*/){rcp=1;}		

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

	double xminus=xepsvle, xplus=1.0-xepsvle, xintminus=0.5,xintplus=0.5;
	double r3=0.5,r4=0.5;
	if(ao_mode==-1)
	{
		r0=0;
		rmax=anzvle*static_cast<double>(dinorm);
		r3=rmax/2.;
		r4=rmax/2;
	}
	if(ao_mode!=-1)
	{
		unordered_map<int,double> *phimtmp=&phi_z_vle[abs(f1)][abs(f2)];

		int ixminus, ixplus,ixintminus,ixintplus;
		vleboundaries( p.Q,  p.th_old, p.x,  ixminus, ixplus,xminus, xplus,ixintminus, ixintplus,xintminus,xintplus);

		double phim2=(*phimtmp).find(ixplus)->second;
		double phim1=(*phimtmp).find(ixminus)->second;

		double phim4=(*phimtmp).find(ixintplus)->second;
		double phim3=(*phimtmp).find(ixintminus)->second;

		double phi1=phmaxvle[abs(f1)][abs(f2)];
	
		double rtmp=anzvle*static_cast<double>(dinorm);

		r3=rtmp*(phim3/phi1);
		r4=rtmp*(phim4/phi1);

		r0=rtmp*(phim1/phi1);
		rmax=rtmp*((phim2-phim1)/phi1);

	}	
	//double r=r0+rand01()*rmax;
	rmax=rmax;//-(r4-r3);
	double r=r0+rand01()*rmax;
	if(r>r3){r=r/*+(r4-r3)*/;}
	unordered_map<int,double> *phinvtmp;

	phinvtmp=&phiinvvle[abs(f1)][abs(f2)];
//	do{
	if(xplus==xminus /*or ixplus==ixminus *//*or (xintminus==xminus and xintplus==xplus)*/){xtmp=0.5;/* cout<<"hi"<<endl;*/}		
	else{
		xtmp=interpolate_integer_inversegrid_1d(phinvtmp, xminus,  xplus,  r);}
//	}while(xtmp>xintminus and xtmp<xintplus);
//	cout<<"xtmp="<<xtmp<<endl;
	return xtmp;
}

void RotateToP(TMDICEparticle p, TMDICEparticle &p1, double phirel)
{
		p1.p0=p1.x*emax;
		p1.px=p1.ktrel*cos(phirel);
		p1.py=p1.ktrel*sin(phirel);
		p1.pz=sqrt(p1.p0*p1.p0-p1.Q*p1.Q-p1.ktrel*p1.ktrel);


		double thtmp=acos(p.pz/sqrt(p.p0*p.p0-p.Q*p.Q)),phtmp=atan2(p.py,p.px);
		double p1ztmp=p1.pz*cos(thtmp)-p1.px*sin(thtmp);
		double p1xtmp=p1.pz*sin(thtmp)+p1.px*cos(thtmp);
		double p1ytmp=p1.py;
		p1.px=p1xtmp*cos(phtmp)-p1ytmp*sin(phtmp);
		p1.py=p1xtmp*sin(phtmp)+p1ytmp*cos(phtmp);
		p1.pz=p1ztmp;
}

void VLEevent::vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2)
{
	//	cout<<"vertex"<<flush<<endl;
	p1.Q_old=p.Q;
	p2.Q_old=p.Q;

	p1.x_old=p.x;
	p2.x_old=p.x;

	p1.typ=typselect(p);	
	p2.typ=typ2[p.typ][p1.typ];
	double kttmp2;
	do{
	p1.x=xselect(p,p1,p2);	
	if(checkdla==0)
	{
		p2.x=1.0-p1.x;//this is now only for selection of z, i.e.: relative fraction
		p1.x=p1.x*p.x;//this is for total energy fraction x
		p2.x=p2.x*p.x;
	}

	if(checkdla==1){p2.x=1.-p1.x;p1.x=1.;p1.typ=p.typ;}
	
	p1.t_old=p.t_old+p1.tf();
	p2.t_old=p.t_old+p2.tf();

	p1.th_old_old=p.th_old;
	p2.th_old_old=p.th_old;
	if(p1.x>p2.x)
	{
		p1.adr=p.adr+"1";
		p2.adr=p.adr+"2";
	}
	else
	{
		p2.adr=p.adr+"1";
		p1.adr=p.adr+"2";		
	}

	p1.Q=evol_select(p1);
	p2.Q=evol_select(p2);

	double z=p1.x/p.x;
	kttmp2=p.Q*p.Q*z*(1.-z)-(1.-z)*p1.Q*p1.Q-z*p2.Q*p2.Q;
	kttmp2=kttmp2-0.25*(pow(p.Q,4.)+pow(p1.Q,4)+pow(p2.Q,4)-2.*pow(p1.Q*p.Q,2)-2.*pow(p.Q*p2.Q,2)-2.*pow(p1.Q*p2.Q,2))/pow(p.x*emax,2);
	kttmp2=kttmp2/(1.-pow(p.Q/(p.x*emax),2.));
	//cout<<p.Q<<" "<<p1.Q<<" "<<p2.Q<<endl;

	p1.ktrel=sqrt(kttmp2);
	p2.ktrel=sqrt(kttmp2);

	p1.kt=p.kt*p1.x;
	p2.kt=p.kt*p2.x;


	p1.phik=p.phik;
	p2.phik=p.phik;	
	p1.th_old=p2.th12();
	p2.th_old=p2.th12();

	double phirel=rand01()*2.*pi;	

	if(rotate==1)
	{
		RotateToP(p, p1, phirel);
		RotateToP(p, p2, phirel+pi);
	}else
	{
		p1.p0=p1.x*emax;
		p1.px=z*p.px+p1.ktrel*cos(phirel);
		p1.py=z*p.py+p1.ktrel*sin(phirel);
		p1.pz=z*p.pz;//sqrt(p1.p0*p1.p0-p1.Q*p1.Q-p1.ktrel*p1.ktrel);

		p2.p0=p2.x*emax;
		p2.px=(1.-z)*p.px+p2.ktrel*cos(phirel+pi);
		p2.py=(1.-z)*p.py+p2.ktrel*sin(phirel+pi);
		p2.pz=(1.-z)*p.pz;//sqrt(p1.p0*p1.p0-p1.Q*p1.Q-p1.ktrel*p1.ktrel);


	}
	p1.kt=sqrt(p1.px*p1.px+p1.py*p1.py);
	p2.kt=sqrt(p2.px*p2.px+p2.py*p2.py);

	p1.phik=atan2(p1.py,p1.px);
	p2.phik=atan2(p2.py,p2.px);
	}while(kttmp2<0 and checkdla!=1);
}