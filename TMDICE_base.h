#include"deps.h"

#ifndef TMDICEbase	
#define TMDICEbase

#include"TMDICE_lib.h"
#include"decor.h"
#include"globalconsts.h"
/**TMDICE_base contains the general structure of the relevant Monte-Carlo algorithms (both VLE and jet medium interactions).
The relevant functions for the creation and evolution of jets are in the abstract class TMDICEbase_event. That class 
contains also the data structures for storing the variables that specify the jet-particles. 
Jet-particles are defined via the class TMDICEparticle inside this file*/

void setparams(map<string,double*> invars2);
/**This structure contains variables to describe an individual jet-particle*/
class TMDICEparticle //maybe throw this away and use only the version that is in TMDICE.h?
{
public:
TMDICEparticle();
double x, t,typ,t_old,kt,phik,Q,th;
double Q_old,x_old,th_old, th_old_old;
bool dump;
string adr;
bool istmp;
bool isfin;
bool islead;

double px,py,pz,p0;
double omega();
double th12();
double kt12(),ktrel,eta(),y(),pt(),phi();
double tf();
double td();

virtual ~TMDICEparticle();
};

//typedef double (*fkt)(double);
typedef double (*fktpart)(TMDICEparticle);
typedef bool (*fktpartbool)(TMDICEparticle);
typedef bool (*fktpartbool3)(TMDICEparticle,TMDICEparticle,TMDICEparticle);
typedef void (*fktpartvoid)(TMDICEparticle&,double);
typedef void (*fktpartvoid0)(TMDICEparticle&);

void rotate4mom(TMDICEparticle &p, double theta, double phi);

int nodumpsize(vector<TMDICEparticle> jet);

/**Abstract class that gives the general structure of the Monte-Carlo algorithm for jet evolution (both for BDIM and VLE) and contains data-structures 
for storing the variables of jet particles (where jet-particles are instances of the TMDICEparticle class).**/
class TMDICEbaseevent
{
public:
TMDICEbaseevent();
vector<TMDICEparticle> gen,genfin,casc;

fktpartbool stopcond, dumpcond;
fktpartbool3 nosplitcond;
fktpartvoid setpval;
fktpart selectionbias;
double evolmin;

void setevolmin(double minq);
void setsetpval(fktpartvoid tmpsetpval);
void setstopcond(fktpartbool tmpstopcond);
void setnosplitcond(fktpartbool3 nosplitcondtmp);
void setdumpcond(fktpartbool tmpdumpcond);
void setselectionbias(fktpart selectionbiastmp);

void setEmax(double ee);
void settmax(double tup);
void settmin(double tmin);

void settyp1(double typ0);
void setkt1(double kt0);
void setx1(double x0); 

virtual double evol_select(TMDICEparticle p)=0;//=0;

virtual double typselect(TMDICEparticle p)=0;

virtual double xselect(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)=0;

virtual void vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2)=0;

virtual void stopparticle(TMDICEparticle p,vector<TMDICEparticle>&newgen,vector<TMDICEparticle>&newgenfin);

virtual void fillgen();

virtual void initgen();

virtual void make_event();

virtual ~TMDICEbaseevent();
};

bool ptrue(TMDICEparticle p);
bool ptrue(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);

bool pfalse(TMDICEparticle p);
bool pfalse(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);

double no_bias(TMDICEparticle p);


#endif
