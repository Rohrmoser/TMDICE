#include "globalconsts.h"
using namespace std;

double gevfm;

double xeps,xmin;
double anz;
double qmin;
double qmax;
double qanz;
double u_min;
double u_max;
double u_anz;

double colorfact,ca,cf;
double tanz=4.*pow(10.,2.);
int time_mode; //time_mode=0...infinite constant medium, 1...static medium of length tmax

double tmin,tmax,taumin,taumax, qhat,qhatgev, emax, alphas, alphabar,nc, ndens, inittyp,W_scat, md,tstar,tstargev,Temp,x0,kt0,typ0; 
double theta0orR,Q2max,xminvle,xepsvle, anzvle, Q2min,ao_mode,R,checkdla, vlemode,simplifythkt,setbias,rotate,p1x,p1y,p1z,bdimmode,vle2mode;
double scatflag1,ktsplitflag1,time_mode1,dlla;

int ktsplitflag, scatflag;

map<string,double> invars;

double pmax=2.;
int dinorm=/*1.137**/pow(10,pmax);

double cqq_u[M], cqg_u[M],cgg_u[M];

map<double,map<double,double>> typ2={{1.0,{{1.0,2.0}}},{2.0,{{1.0,-1.0},{-1.0,1.0},{2.0,2.0}} },{-1.0,{{-1.0,2.0}}} };

void readTMDICEparameters(string snom)
{
	disclaimer();
	string comment;
	ifstream s;
	cout<<endl<<"This is the name of the input file: "<<snom<<endl;
	s.open(snom.c_str() );
	cout<<endl<<"Input parameters:"<<endl;
	while(s.eof()==false )
	{
		string vnam;
		double vval;
		s>>vnam>>vval;
		getline(s,comment);
		if(vnam!="")
			{
				invars[vnam]=vval;
				cout<<vnam<<" "<<invars[vnam]<<comment<<endl;
			}
	}
	s.close();
	string errmess="The following information is necessary: 1.) nc and ndens, qhat and emax 2.)Either taumin and taumax or tmin and tmax 3.) alphas or alphabar";
	auto fin=invars.end();
	if(invars.find("nc")==fin or invars.find("ndens")==fin or 
	invars.find("qhat")==fin or invars.find("emax")==fin ){cout<<errmess<<endl; exit(1);}
	if((invars.find("taumin")==fin or invars.find("taumax")==fin) and (invars.find("tmin")==fin or invars.find("tmax")==fin)){cout<<errmess<<endl;exit(2);}
	if(invars.find("alphas")==fin and invars.find("alphabar")==fin){cout<<errmess<<endl;exit(3);}
}

void readTMDICEparameters(map<string,double>iv1)
{
	disclaimer();
	
	invars=iv1;
	string errmess="The following information is necessary: 1.) nc and ndens, qhat and emax 2.)Either taumin and taumax or tmin and tmax 3.) alphas or alphabar";
	auto fin=invars.end();
	if(invars.find("nc")==fin or invars.find("ndens")==fin or 
	invars.find("qhat")==fin or invars.find("emax")==fin ){cout<<errmess<<endl; exit(1);}
	if((invars.find("taumin")==fin or invars.find("taumax")==fin) and (invars.find("tmin")==fin or invars.find("tmax")==fin)){cout<<errmess<<endl;exit(2);}
	if(invars.find("alphas")==fin and invars.find("alphabar")==fin){cout<<errmess<<endl;exit(3);}
}

map<double,double> up_lim_phi_time={{0,1},{1,sinh(pi)/(cosh(pi)-1.)},{2,1.5}};

