#include"deps.h"
#include"TMDICE_lib.h"
#include"decor.h"
#include "globalconsts.h"
#include "functions_lib.h"
#include"TMDICE_base.h"

using namespace std;

void setparams(map<string,double*> invars2)
{	
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
double TMDICEparticle::omega()
{
	double x_old=old_part->x;

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



/**@brief branching angle

 of the current particle with respect to its parent.
*/
double TMDICEparticle::th12()
{
	double x_old=old_part->x;
	double Q_old=old_part->Q;

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

/**@brief transverse momentum 

of the current jet particle with respect to its parent.
*/
double TMDICEparticle:: kt12()
{
	double kttmp;
	if(checkdla==0)
	{
		kttmp=ktrel;
	}
	if(checkdla==1 or simplifythkt==1)
	{
		kttmp=omega()*th_old;
	}
	return kttmp;
}

double TMDICEparticle:: tf()
{
	double x_old=old_part->x;
	double Q_old=old_part->Q;
	return (2.*emax*x_old/(Q_old*Q_old))/gevfm;
}

double TMDICEparticle:: td()
{
	double thold=th_old;
	return pow(12./(qhatgev*thold*thold),1./3.)/gevfm;
}

void TMDICEbaseevent::stopparticle(TMDICEparticle *p,vector<TMDICEparticle*>& cur,vector<TMDICEparticle*>& fin)
{
	if(dumpcond(*p)==false)
	{
		p->dump=false;
		if(stopcond(*p)){ p->isfin=true;fin.push_back(p);}
		else{cur.push_back(p);}

	}
	else{p->dump=true;setpval(*p,evolmin); p->isfin=true;if(dumpthedump ){p=NULL;}else{fin.push_back(p);}}
}

void  TMDICEbaseevent::fillgen(vector<TMDICEparticle*> &cur,vector<TMDICEparticle*>&fin)
{
	if(cur.size()>0)
	{
		vector<TMDICEparticle*>genptmp=cur;
		cur.clear();
		for( int j=0; j<genptmp.size(); j++)
		{
			TMDICEparticle *p=genptmp.at(j);

			TMDICEparticle p1,p2;
			
			vertex(*p,p1,p2);
			p->insert(p1);
			p->insert(p2);
			if(nosplitcond(*p, *(p->left),*(p->right) ) )
			{
				if(filltemps)
				{
					p->left->isfin=true;
					fin.push_back(p->left);
					p->right->istmp=true;
				}
				else
				{
					p->isfin=true;
					fin.push_back(p);

					delete p->left;
					delete p->right;
					p->left=NULL;
					p->right=NULL;
				}
			}
			else
			{				
				stopparticle(p->left,cur,fin);
				stopparticle(p->right,cur,fin);	
			}		
		}
	}
}

void TMDICEbaseevent::initgen(TMDICEparticle &c,vector<TMDICEparticle*> &cur,vector<TMDICEparticle*>&fin)
{
	cur.clear();
	double varpi=pi;

	fourmom tmpmom(0.,0.,0.,0.);
	TMDICEparticle p0(x0,tmin,typ0,kt0,0.,0.,min(varpi,R),false,"1",false,false, true,tmpmom);
	TMDICEparticle p00(x0,tmin,typ0,kt0,0.,Q2max,min(varpi,R),false,"0",false,false, true,tmpmom);
	c=p00;
	c.insert(p0);

	TMDICEparticle * p0tmp=c.left;
	double qtmp=evol_select(*p0tmp);
	setpval(*p0tmp,qtmp);
	p0tmp->p.p0=p0.x*emax;
	double pp=sqrt(p0tmp->p.p0*p0tmp->p.p0-p0tmp->Q*p0tmp->Q)/p0tmp->p.p0;
	
	p0tmp->p.px=p1x*pp;
	p0tmp->p.py=p1y*pp;
	p0tmp->p.pz=p1z*pp;

	stopparticle(c.left,cur,fin);
}


void TMDICEbaseevent::make_event(TMDICEparticle &c,vector<TMDICEparticle*> &cur,vector<TMDICEparticle*>&fin)
{
	initgen(c,cur,fin);
	while(cur.size()>0)
	{
		fillgen(cur,fin);
	}
}

bool fincond(TMDICEparticle *p)
{
    return (p->isfin)==true;
}

bool tempcond(TMDICEparticle *p)
{
    return (p->istmp)==true;
}
bool ptrue(TMDICEparticle p)
{
	return true;
}

bool ptrue(TMDICEparticle *p)
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

void determineadress(TMDICEparticle& p1, TMDICEparticle& p2)
{
	TMDICEparticle *p=p1.old_part;
	if(p1.x>p2.x)
	{
		p1.adr=p->adr+"1";
		p2.adr=p->adr+"2";
		p2.islead=0;
	}
	else
	{
		p2.adr=p->adr+"1";
		p1.adr=p->adr+"2";
		p1.islead=0;		
	}
}
