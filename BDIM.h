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
//void setBDIM();

void setkernelsbdim();

void obtaincdfsbdim();

//dumpcond
bool ltxmin(TMDICEparticle p);
//stopcond
bool stop_with_time(TMDICEparticle p);
//setpval
void set_to_t(TMDICEparticle &p, double tt);

/**BDIMevent is a class that has as its output-data vectors of particles (the structure TMDICEparticle).
These vectors are genfin (for the final particles at time tmax) and casc (for all the particles; containing intermediate particles as well.
There also is an auxiliary vector gen, which describes the set of unevolved/intermediate particles, which should be empty at the end of the evolution. */
class BDIMevent : public TMDICEbaseevent
{
public:

BDIMevent();
void scat(TMDICEparticle p, TMDICEparticle &p1);

string interact_select(TMDICEparticle p);
double evol_select(TMDICEparticle p);

double typselect(TMDICEparticle p);

double xselect(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);

void vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2);

void fillgen(vector<TMDICEparticle*> &cur,vector<TMDICEparticle*> &fin);
};


#endif
