#include"deps.h"
#include"TMDICE_lib.h"
#include"decor.h"
#include "globalconsts.h"
#include "functions_lib.h"
#include"TMDICE_base.h"


void setparams(map<string,double*> invars2)
{	
	//double scatflag1,ktsplitflag1,time_mode1;
	for(map<string,double*>::iterator it=invars2.begin();it!=invars2.end();it++)
	{
		if(invars.find(it->first)!=invars.end())
		{
			cout<<invars.find(it->first)->first<<"=";
			double *tmp=it->second;
			*tmp=invars[it->first];
			cout<<*tmp<<endl;
		}
	}
	scatflag=static_cast<int>(scatflag1);
	ktsplitflag=static_cast<int>(ktsplitflag1);
	time_mode=static_cast<int>(time_mode1);
}
//TMDICEparticle class:
TMDICEparticle::TMDICEparticle()
{
	
}

double TMDICEparticle::omega()
{
	double z=x/x_old;
	if(checkdla==1)
	{
		return z*emax;
	}
	else
	{
		return z*(1.0-z)*x_old*emax;
	}	
}

double TMDICEparticle::eta()
{
	double p=sqrt(px*px+py*py+pz*pz);
	return 0.5*log((p+pz)/(p-pz));
}

double TMDICEparticle::y()
{
	return 0.5*log((p0+pz)/(p0-pz));
}

double TMDICEparticle::phi()
{
	return atan2(py,px);
}

double TMDICEparticle::pt()
{
	return sqrt(px*px+py*py);
}


double TMDICEparticle::th12()
{
	double z=x/x_old;

	double thtmp=1.-(Q_old*Q_old)/(2.*x_old*emax*omega());
	if(checkdla==0)
	{
		thtmp=acos(thtmp);
	}
	if(checkdla==1)
	{
		thtmp=sqrt( (Q_old*Q_old)/(emax*omega()));
	}
	if(simplifythkt==1)
	{
		thtmp=(Q_old/(x_old*emax))/sqrt(z*(1.-z));
	}
	return thtmp;
}

double TMDICEparticle:: kt12()
{
	double kttmp;
	if(checkdla==0)
	{
		kttmp=ktrel;
	}
	if(checkdla==1 or simplifythkt==1)
	{
		kttmp=omega()*th12();
	}
	return kttmp;
}

double TMDICEparticle:: tf()
{
	return (2.*emax*x_old/(Q_old*Q_old))/gevfm;
}

double TMDICEparticle:: td()
{
	double thold=th_old;
	return pow(12./(qhatgev*thold*thold),1./3.)/gevfm;
}

TMDICEparticle::~TMDICEparticle()
{

}

//TMDICEbaseevent class:
TMDICEbaseevent::TMDICEbaseevent()
{
	
}

void TMDICEbaseevent::setEmax(double ee)
{
	emax=ee;
}
void TMDICEbaseevent::settmax(double tup)
{
	tmax=tup;
}

void TMDICEbaseevent::settmin(double tdown)
{
	tmin=tdown;
}

void TMDICEbaseevent::setx1(double x00)
{
	x0=x00;
}

void TMDICEbaseevent::setkt1(double kt00)
{
	kt0=kt00;
}

void TMDICEbaseevent::settyp1(double typ00)
{
	typ0=typ00;
}

void rotate4mom(TMDICEparticle &p, double theta, double phi)
{
	double pabs=sqrt(p.px*p.px+p.py*p.py+p.pz*p.pz);
	p.px=pabs*sin(theta)*cos(phi);
	p.px=pabs*sin(theta)*sin(phi);
	p.px=pabs*cos(theta);
}


void TMDICEbaseevent::stopparticle(TMDICEparticle p, vector<TMDICEparticle>&newgen,vector<TMDICEparticle>&newgenfin)
{
	if(dumpcond(p)==false)
	{
		p.dump=false;
		if(stopcond(p)){ newgenfin.push_back(p);}
		else{newgen.push_back(p);}

	}
	else{p.dump=true;setpval(p,evolmin); newgenfin.push_back(p);}
}

void  TMDICEbaseevent::fillgen()
{
	vector<TMDICEparticle> newgen,newgenfin,newgenfin2;
	if(gen.size()>0)
	{
		for( int j=0; j<gen.size(); j++)
		{
			TMDICEparticle p=gen.at(j);
			

			TMDICEparticle p1,p2;
			
			vertex(p,p1,p2);
			if(nosplitcond(p, p1,p2))
			{
				newgenfin2.push_back(p);
			}
			else
			{				
				stopparticle(p1,newgen,newgenfin);
				stopparticle(p2,newgen,newgenfin);	
			}		
		}
	}
	gen=newgen;
	genfin.insert(genfin.end(),newgenfin.begin(),newgenfin.end());
	genfin.insert(genfin.end(),newgenfin2.begin(),newgenfin2.end());
	casc.insert(casc.end(),gen.begin(),gen.end());
	casc.insert(casc.end(),newgenfin.begin(),newgenfin.end());
}


void TMDICEbaseevent::initgen()
{
	vector<TMDICEparticle>newgen;
	TMDICEparticle p0;
	p0.adr="1";
	p0.x=x0;
	p0.x_old=x0;

	p0.typ=typ0;
	p0.t_old=tmin;
	p0.t=tmin;
	//p0.Q=Q2max;
	p0.Q=0;
	double varpi=pi;
	p0.th_old=min(varpi,R);
	p0.th_old_old=p0.th_old;
	p0.Q_old=Q2max;
	p0.ktrel=0;
//	cout<<"hi"<<flush<<endl;
	double qtmp=evol_select(p0);
//	cout<<"hi again"<<endl;
	setpval(p0,qtmp);
	//missing: conversion Q to t
	p0.kt=kt0;
	p0.phik=0;
	p0.islead=1;

	/*p0.p0=p0.x*emax;
	p0.pz=p0.kt;
	p0.py=0.;
	p0.px=sqrt(p0.p0*p0.p0-p0.kt*p0.kt-p0.Q*p0.Q);*/


	p0.p0=p0.x*emax;

	double pp=sqrt(p0.p0*p0.p0-p0.Q*p0.Q)/p0.p0;
	p0.px=p1x*pp;
	p0.py=p1y*pp;
	p0.pz=p1z*pp;

	p0.isfin=false;
//	cout<<p0.p0<<" "<<p0.px<<" "<<p0.py<<" "<<p0.pz<<endl;


	stopparticle(p0,newgen,genfin);
	gen=newgen;
	casc.insert(casc.end(),gen.begin(),gen.end());
	casc.insert(casc.end(),genfin.begin(),genfin.end());
//	cout<<"First particle selected"<<std::flush<<endl;
}

void TMDICEbaseevent::make_event()
{
	gen.clear();
	genfin.clear();
	casc.clear();
	initgen();
	while(gen.size()>0)
	{
		fillgen();
	}
}
void TMDICEbaseevent::setevolmin(double minq)
{
	evolmin=minq;
}

void TMDICEbaseevent::setsetpval(fktpartvoid tmpsetpval)
{
	setpval=tmpsetpval;
}

void TMDICEbaseevent::setstopcond(fktpartbool tmpstopcond)
{
	stopcond=tmpstopcond;
}

void TMDICEbaseevent::setnosplitcond(fktpartbool3 nosplitcondtmp)
{
	nosplitcond=nosplitcondtmp;
}

void TMDICEbaseevent::setdumpcond(fktpartbool tmpdumpcond)
{
	dumpcond=tmpdumpcond;
}

void TMDICEbaseevent::setselectionbias(fktpart selectionbiastmp)
{
	selectionbias=selectionbiastmp;
}

TMDICEbaseevent::~TMDICEbaseevent()
{
	
}

bool ptrue(TMDICEparticle p)
{
	return true;
}

bool ptrue(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	return true;
}

bool pfalse(TMDICEparticle p)
{
	return false;
}

bool pfalse(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	return false;
}

double no_bias(TMDICEparticle p)
{
	return -0.5;
}

int nodumpsize(vector<TMDICEparticle> jet)
{
	int j=0;
	for(int i=0;i<jet.size();i++)
	{
		if(jet.at(i).dump==false)
		{
			j++;
		}
	}
	return j;
}