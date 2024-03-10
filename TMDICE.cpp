#include"deps.h"

#include"TMDICE_lib.h"
#include"decor.h"
#include"globalconsts.h"
#include "TMDICE_base.h"
#include "BDIM.h"
#include "vle.h"
#include "TMDICE.h"
#include "functions_lib.h"


using namespace std;

void setdefaultvalues()
{
    gevfm=1.0/0.19732;	
	xeps=1.*pow(10,-4);
	xmin=1.*pow(10,-4);
	anz=2*pow(10,7);

	qmin=.1;
	qmax=40.;	
	qanz=2*pow(10,4);	
	u_min=0.;	
	u_max=pi;
	u_anz=1*pow(10.,7);

	kt0=0.;
	x0=1.;
	typ0=2;
	N_F=3.;
	scatflag=0;
	ktsplitflag=1;
	time_mode=0;
	rotate=1;

    bdimmode=1;
    vlemode=0;
    vle2mode=1;

    //now parameters for vle
	xepsvle=1.*pow(10,-4);
	xminvle=1.*pow(10,-4);
	anzvle=3.*pow(10,6);
	ao_mode=1;
	R=pi;
	checkdla=0;
	setbias=0;
	simplifythkt=0;
}

map<string,double*> invars2=
{{"xeps",&xeps},{"xmin",&xmin},{"qmin",&qmin},{"kt1",&kt0},{"x1",&x0},{"typ1",&typ0},{"nc",&nc},
{"alphabar",&alphabar},{"alphas",&alphas},{"ndens",&ndens},{"T",&Temp},{"qhat",&qhat},{"emax",&emax},
{"md",&md},{"taumin",&taumin},{"taumax",&taumax},{"tmin",&tmin},{"tmax",&tmax},{"scat",&scatflag1},{"ktsplit",&ktsplitflag1}
,{"timemode",&time_mode1},{"rotate",&rotate},{"p1x",&p1x},{"p1y",&p1y},{"p1z",&p1z},
//now for VLE:
{"xepsvle",&xepsvle},{"xminvle",&xminvle},{"theta0orR",&theta0orR},{"R",&R},{"Q2max",&Q2max},{"Q2min",&Q2min},
{"ao_mode",&ao_mode},{"ktmin",&ktmin},{"checkdla",&checkdla},{"simplifythkt",&simplifythkt},{"setbias",&setbias},{"p1x",&p1x}
,{"p1y",&p1y},{"p1z",&p1z},{"rotate",&rotate},
{"vlemode",&vlemode},{"vle2mode",&vle2mode},{"bdimmode",&bdimmode},{"dlla",&dlla}
};

void derivedvalues()
{
    ca=nc;
	cf=(nc*nc-1.)/(2.*nc);

	colorfact=pow(cf/ca,1.0);

	if(invars.find("alphabar")==invars.end())
	{
		alphabar=invars["alphas"]/pi;
		alphas=invars["alphas"];
	}
	else
	{
		alphabar=invars["alphabar"];
		alphas=alphabar*pi;
	}
	qhatgev=qhat/gevfm;//inGeV^3
	tstar=alphabar*sqrt(qhatgev/emax);//inGeV
	tstargev=(1.0/tstar);
	tstar=tstargev/gevfm;//infm

	if(invars.find("taumin")!=invars.end() and invars.find("taumax")!=invars.end() ) 
	{
		tmin=taumin*tstar;
		tmax=taumax*tstar;
	}
	if(invars.find("tmin")!=invars.end() and invars.find("tmax")!=invars.end() ) 
	{
		taumin=tmin/tstar;
		taumax=tmax/tstar;
	}
}

void setkernels()
{
    if(bdimmode==1)
    {
        setkernelsbdim();
    }

}

void obtaincdfs()
{
    if(bdimmode==1)
    {
        obtaincdfsbdim();
    }

    if(vlemode==1)
    {
        cout<<"obtain CDFs for VLE:"<<endl;
        obtaincdfsvle();
    }
}

void setTMDICE()
{
	disclaimer();
    setdefaultvalues();
	setparams(invars2);
    derivedvalues();
    setkernels();
	cout<<endl<<"Calculating partition functions: This may take a while."<<endl;
    obtaincdfs();
	srand (time(NULL));
	cout<<"Initialization finished"<<endl;
}



double TMDICEevent::evol_select(TMDICEparticle p){return 0.;}
double TMDICEevent::typselect(TMDICEparticle p){return 0.;}
double TMDICEevent::xselect(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2){return 0.;}
void TMDICEevent::vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2){}


void TMDICEevent::make_event()
{

    VLEevent vle1(stop_non_resolved_but_resolvable,stop_non_resolved_but_resolvable);
    VLEevent vle2(stop_with_ktmin,stop_with_ktmin);
    BDIMevent bdim;
    bdim.dumpthedump=dumpthedump;
    vle1.dumpthedump=dumpthedump;
    vle2.dumpthedump=dumpthedump;
    bdim.cascfill=false;
    vle1.cascfill=false;
    vle2.cascfill=false;

    vector<TMDICEparticle*> g1,g2;
    g1.clear();
    g2.clear();
    if(vlemode)
    {
        vle1.make_event(root,g1,g2);
    }
    else
    {
        bdim.initgen(root,g1,g2);
    }
    
    if(bdimmode)
    {
        
        vector<TMDICEparticle*> g3=g2;
        g3.clear();
        g3.insert(g3.end(),g2.begin(),g2.end());

        g3.insert(g3.end(),g1.begin(),g1.end());

        g1.clear();
        g2.clear();
        for(int i=0; i<g3.size();i++)
        {
            g3.at(i)->isfin=false;
            bdim.stopparticle(g3.at(i),g1,g2);
        }
        while(!g1.empty())
        {
            bdim.fillgen(g1,g2);
        }
    }

    if(vlemode and vle2mode)
    {
        g1.clear();
        for(int i=0; i<g2.size();i++)
        {
            g2.at(i)->isfin=false;
            if(!vle2.dumpcond(*g2.at(i))and g2.at(i)->Q>Q2min){g1.push_back(g2.at(i));}
        }
        g2.clear();
        while(!g1.empty())
        {
            vle2.fillgen(g1,g2);
        }
    }
    genp=g1;
    genfinp.clear();
    vector<TMDICEparticle*> all=getcasc();
    for(int i=0;i<all.size();i++)
    {
        if(all.at(i)->isfin)
        {
            genfinp.push_back(all.at(i));
        }
    } 
}

