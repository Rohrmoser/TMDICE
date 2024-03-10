#include "deps.h"

#include "TMDICE_lib.h"
#include "decor.h"
#include "globalconsts.h"
#include "TMDICE_base.h"
#include "BDIM.h"
#include "vle.h"
#ifndef TMDICE	
#define TMDICE

void setTMDICE();
void derivedvalues();

class TMDICEevent : public TMDICEbaseevent
{
public:
TMDICEevent(){};

double evol_select(TMDICEparticle p);

double typselect(TMDICEparticle p);

double xselect(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);

void vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2);

void make_event();
};
#endif
