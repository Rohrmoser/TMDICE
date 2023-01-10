#include"deps.h"

#include"TMDICE_lib.h"
#include"decor.h"
#include"globalconsts.h"
#include "TMDICE_base.h"
#include "BDIM.h"
#include "vle.h"
#include "TMDICE.h"

void setTMDICE()
{
    vlemode=0;
    if(invars.find("vlemode")!=invars.end()){vlemode=invars["vlemode"];}
    if(vlemode==1)
    {
        setVLE();
    }
    setBDIM();
}



double TMDICEevent::evol_select(TMDICEparticle p){return 0.;}
double TMDICEevent::typselect(TMDICEparticle p){return 0.;}
double TMDICEevent::xselect(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2){return 0.;}
void TMDICEevent::vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2){}


void TMDICEevent::make_event()
{
    if(vlemode==0)
    {
        BDIMevent jetbdim;
        jetbdim.make_event();
        casc=jetbdim.casc;
        gen=jetbdim.gen;
        genfin=jetbdim.genfin;
    }
    if(vlemode==1)
    {
        BDIMevent bdim;
        VLEevent vle1(stop_non_resolved_but_resolvable,stop_non_resolved_but_resolvable);
        VLEevent vle2(stop_with_ktmin,stop_with_ktmin);
        vle1.make_event();
        bdim.casc=vle1.casc;
        bdim.gen.clear();
        bdim.genfin.clear();
        //cout<<"obtained VLE1"<<endl;
        //check if particles can still evolve:
        for(int i=0; i<vle1.genfin.size();i++)
        {
            TMDICEparticle p=vle1.genfin.at(i);
            bdim.stopparticle(p,bdim.gen, bdim.genfin);
        }
        while(bdim.gen.size()>0)
        {
            bdim.fillgen();
        }
       // cout<<"obtained BDIM"<<endl;
        vle2.casc=bdim.casc;
        vle2.gen.clear();
        vle2.genfin.clear();
        //check if particles can still evolve:
        for(int i=0; i<bdim.genfin.size();i++)
        {
            TMDICEparticle p=bdim.genfin.at(i);
            vle2.stopparticle(p,vle2.gen, vle2.genfin);
        }
        while(vle2.gen.size()>0)
        {
            vle2.fillgen();
        }    
     //   cout<<"obtained VLE2"<<endl;    
        casc=vle2.casc;
        gen=vle2.gen;
        genfin=vle2.genfin;
    }
}

