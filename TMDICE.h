#include"deps.h"

#ifndef TMDICE	
#define TMDICE

#include"TMDICE_lib.h"

void readTMDICEparameters(string snom);
void readTMDICEparameters(map<string,double>iv1);


void setTMDICE();

struct TMDICEparticle 
{
double x, t,typ,t_old,kt,phik;
bool dump;
};

class TMDICEevent
{
public:

TMDICEevent();

void setEmax(double ee);
void settmax(double tup);
void settmin(double tmin);

void settyp1(double typ0);
void setkt1(double kt0);
void setx1(double x0);

vector<TMDICEparticle> gen,genfin,casc;

void fillgen();

void initgen();

double tselect(double f,double x, double told);


double typselect(double f1);

double xselect(double f1, double f2);

string interact_select(double f,double x);


void make_event();

virtual ~TMDICEevent();
};


#endif
