#include"deps.h"
#include"TMDICE_lib.h"
#include"decor.h"
#include"globalconsts.h"
#include"TMDICE_base.h"

#ifndef VLE	
#define VLE
extern double ktmin;
/**
gives the parameters for the TMDICE program to run properly to the TMDICE program and calculates the necessary cumulative distribution functions.
*/
void setVLE();

bool stop_with_ktmin(TMDICEparticle p);
bool stop_with_ktmin(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);
void settoq(TMDICEparticle &p,double q);
bool ltxminvle(TMDICEparticle p);
bool stop_with_qmin(TMDICEparticle p);
bool stop_when_tf_exceeds_td(TMDICEparticle p);
bool stop_when_td_exceeds_tmax(TMDICEparticle p);
bool stop_when_tf_exceeds_tmax(TMDICEparticle p);
bool stop_when_tf_exceeds_tmax_and_ktmin(TMDICEparticle p);
bool stop_when_td_exceeds_tmax_and_ktmin(TMDICEparticle p);
bool stop_non_resolved_but_resolvable(TMDICEparticle p);

bool stop_when_tf_exceeds_td(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);
bool stop_when_td_exceeds_tmax(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);
bool stop_when_td_exceeds_tmax_and_ktmin(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);
bool stop_when_tf_exceeds_tmax(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);
bool stop_when_tf_exceeds_tmax_and_ktmin(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);
bool stop_non_resolved_but_resolvable(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);
double bias_ktmin(TMDICEparticle p);
double bias_tftd(TMDICEparticle p);
double bias_tftd_ktmin(TMDICEparticle p);
double bias_tdtl(TMDICEparticle p);

class VLEevent : public TMDICEbaseevent
{
public:

VLEevent();
VLEevent( fktpartbool stopcondtmp, fktpartbool3 nosplitcondtmp=pfalse,fktpart selectionbiastmp=no_bias,fktpartvoid pvaltmp=settoq, fktpartbool dumpcondtmp=ltxminvle,double minq=Q2min);

double evol_select(TMDICEparticle p);

double typselect(TMDICEparticle p);

double xselect(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);

void vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2);

void vleboundaries(double q, double thold,double x, int & ixminus, int& ixplus, double & xminus, double & xplus,
 int & ixintminus, int& ixintplus, double & xintminus, double & xintplus);
};

#endif