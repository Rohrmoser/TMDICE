#include"deps.h"

#ifndef BDIM	
#define BDIM

#include"TMDICE_lib.h"
#include"decor.h"
#include"globalconsts.h"
#include"TMDICE_base.h"


/**
gives the parameters for the TMDICE program to run properly to the TMDICE program and calculates the necessary cumulative distribution functions.
*/
void setBDIM();

//dumpcond
bool ltxmin(TMDICEparticle p);
//stopcond
bool stop_with_time(TMDICEparticle p);
//setpval
void set_to_t(TMDICEparticle &p, double tt);

/**TMDICEevent is a class that has as its output-data vectors of particles (the structure TMDICEparticle).
These vectors are genfin (for the final particles at time tmax) and casc (for all the particles; containing intermediate particles as well.
There also is an auxiliary vector gen, which describes the set of unevolved/intermediate particles, which should be empty at the end of the evolution. */
class BDIMevent : public TMDICEbaseevent
{
public:

//TMDICEevent();

//void setEmax(double ee);
//void settmax(double tup);
//void settmin(double tmin);

//void settyp1(double typ0);
//void setkt1(double kt0);
//void setx1(double x0);

//vector<TMDICEparticle> gen,genfin,casc;

//void med_vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2);

//need a neat constructor
BDIMevent();
//BDIMevent( fktpartbool stopcondtmp=stop_with_time, fktpartvoid pvaltmp=set_to_t, fktpartbool dumpcondtmp=ltxmin,double minq=tmax);

void scat(TMDICEparticle p, TMDICEparticle &p1);

string interact_select(TMDICEparticle p);
double evol_select(TMDICEparticle p);

double typselect(TMDICEparticle p);

double xselect(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);

void vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2);


void fillgen();
/*void fillgen2();

void initgen();
void initgen2();


double tselect(double f,double x, double told);


double typselect(double f1,double t,double x);

double xselect(double f1, double f2);
double xselect(double f1, double f2, double t);

string interact_select(double f,double x,double t);


void make_event();
void make_event2();


virtual ~TMDICEevent();*/
};


#endif
