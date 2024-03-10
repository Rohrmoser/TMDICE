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
/**This class contains variables to describe an individual jet-particle*/
class TMDICEparticle 
{
public:
double x, t,typ,kt,phik,Q,th,th_old;
double kt_sp;
bool dump;
/**@brief individual adress for the particle within the cascade.

The initial particle is 0. For every splitting a number, either "1" or "2" is added, 
identifying the particle as the first or second particle emitted in the splitting.
*/
string adr;
bool istmp;
bool isfin;
bool islead;

fourmom p;
double omega();
double th12();
double kt12(),ktrel;
double eta(){return p.eta();};
double y(){return p.y();};
double pt(){return p.pt();};
double phi(){return p.phi();};
double tf();
double td();

TMDICEparticle* old_part;

TMDICEparticle *left;
TMDICEparticle *right;

TMDICEparticle()
{
    old_part=NULL;
    left=NULL;
    right=NULL;
}

TMDICEparticle(const TMDICEparticle& org)
{
	old_part=org.old_part;
    x=org.x;
    t=org.t;
    typ=org.typ;
    kt=org.kt;
    phik=org.phik;
    Q=org.Q;
    th_old=org.th_old;
    dump=org.dump;
    adr=org.adr;
    istmp=org.istmp;
    isfin=org.isfin;
    islead=org.islead;
    p=org.p;
    ktrel=org.ktrel;

    if(org.left!=NULL)
    {
        left=new TMDICEparticle(*(org.left));
    }
    else
    {
        left=NULL;
    }

    if(org.right!=NULL)
    {
        right=new TMDICEparticle(*(org.right));
    }
    else
    {
        right=NULL;
    }
    reassign();
}

TMDICEparticle(double x2,double t2,double typ2,double kt2,double phik2,double Q2,double th_old2,
bool dump2,string adr2,bool istmp2,bool isfin2, bool islead2,fourmom p2)
{
    old_part=NULL;
    left=NULL;
    right=NULL;
    x=x2;
    t=t2;
    typ=typ2;
    kt=kt2;
    phik=phik2;
    Q=Q2;
    th_old=th_old2;
    dump=dump2;
    adr=adr2;
    istmp=istmp2;
    isfin=isfin2;
    islead=islead2;
    p=p2;
}

TMDICEparticle& operator=(const TMDICEparticle& rhs )
{
    if(this!=&rhs)
    {
        old_part=rhs.old_part;
        x=rhs.x;
        t=rhs.t;
        typ=rhs.typ;
        kt=rhs.kt;
        phik=rhs.phik;
        Q=rhs.Q;
        th_old=rhs.th_old;
        dump=rhs.dump;
        adr=rhs.adr;
        istmp=rhs.istmp;
        isfin=rhs.isfin;
        islead=rhs.islead;
        p=rhs.p;
        ktrel=rhs.ktrel;

        if(rhs.left!=NULL)
        {
            left=new TMDICEparticle(*(rhs.left));
        }
        else
        {
            left=NULL;
        }

        if(rhs.right!=NULL)
        {
            right=new TMDICEparticle(*(rhs.right));
        }
        else
        {
            right=NULL;
        }
        this->reassign();
    }
    return *this;
}

void reassign(TMDICEparticle* here)
{
    if(here->left!=NULL)
    {
        here->left->old_part=here;
        reassign(here->left);
    }
    if(here->right!=NULL)
    {
        here->right->old_part=here;
        reassign(here->right);
    }
}
void reassign()
{
    reassign(this);
}

void insert(const TMDICEparticle& rhs)
{
    inserthere(this,rhs);
}

void inserthere(TMDICEparticle* here,const TMDICEparticle& rhs )
{
    if(this!=&rhs and this->left==NULL)
    {
        this->left=new TMDICEparticle(rhs);
		this->left->old_part=here;

    }
    else
    {
        if(this!=&rhs and this->right==NULL)
        {
            this->right=new TMDICEparticle(rhs);
			this->right->old_part=here;

        }
        else
        {
            inserthere(this->left,rhs);
        }
    }
}

TMDICEparticle* getroot()
{
    TMDICEparticle* root=this;
    while (root->old_part!=NULL){root=root->old_part;}
    return root;
}

vector<TMDICEparticle*> getnextgen(vector<TMDICEparticle*>gen)
{
    vector<TMDICEparticle*>nextgen;
    nextgen.clear();
    for(int i=0; i<gen.size();i++)
    {
        TMDICEparticle* act=gen.at(i);
        if(act->left!=NULL){nextgen.push_back(act->left);}
        if(act->right!=NULL){nextgen.push_back(act->right);}
    }
    return nextgen;
}

vector<TMDICEparticle*> find(bool (*cond)(TMDICEparticle*))
{
    vector<TMDICEparticle*> all=getall(),res;
    res.clear();
    for(int i=0; i<all.size();i++)
    {
        TMDICEparticle* act=all.at(i);
        if(cond(act))
        {
            res.push_back(act);
        }
    }
    return res;
}

vector<TMDICEparticle*> getall()
{
    TMDICEparticle* root=this->getroot();
    vector<TMDICEparticle*>res,res2;
    res.clear();
    res.push_back(root);
    res2=res;
    while(!res2.empty())
    {
        res2=getnextgen(res2);

        res.insert(res.end(),res2.begin(),res2.end());
    }
    return res;
}


virtual ~TMDICEparticle()
{
    delete left;
    delete right;
}

};

typedef double (*fktpart)(TMDICEparticle);
typedef bool (*fktpartbool)(TMDICEparticle);
typedef bool (*fktpartbool3)(TMDICEparticle,TMDICEparticle,TMDICEparticle);
typedef void (*fktpartvoid)(TMDICEparticle&,double);
typedef void (*fktpartvoid0)(TMDICEparticle&);

bool fincond(TMDICEparticle *p);
bool tempcond(TMDICEparticle *p);


bool ptrue(TMDICEparticle p);
bool ptrue(TMDICEparticle *p);

bool ptrue(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);

bool pfalse(TMDICEparticle p);
bool pfalse(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2);

/**Abstract class that gives the general structure of the Monte-Carlo algorithm for jet evolution (both for BDIM and VLE) and contains data-structures 
for storing the variables of jet particles (where jet-particles are instances of the TMDICEparticle class).**/
class TMDICEbaseevent
{
public:
bool dumpthedump=false;
bool cascfill=true;
TMDICEparticle root;
vector<TMDICEparticle*> genp, genfinp;
TMDICEbaseevent(){};

vector<TMDICEparticle> gen,genfin,casc,temps;
bool filltemps;

fktpartbool stopcond, dumpcond;
fktpartbool3 nosplitcond;
fktpartvoid setpval;
double evolmin;

void setevolmin(double minq){evolmin=minq;};
void setsetpval(fktpartvoid tmpsetpval){setpval=tmpsetpval;};
void setstopcond(fktpartbool tmpstopcond){stopcond=tmpstopcond;};
void setnosplitcond(fktpartbool3 nosplitcondtmp){nosplitcond=nosplitcondtmp;};
void setdumpcond(fktpartbool tmpdumpcond){dumpcond=tmpdumpcond;};

void setEmax(double ee){emax=ee;};
void settmax(double tup){tmax=tup;};
void settmin(double tdown){tmin=tdown;};
void setx1(double x00){x0=x00;};
void setkt1(double kt00){kt0=kt00;};
void settyp1(double typ00){typ0=typ00;};

virtual double evol_select(TMDICEparticle p)=0;

virtual double typselect(TMDICEparticle p)=0;

virtual double xselect(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)=0;

virtual void vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2)=0;

virtual void stopparticle(TMDICEparticle *p,vector<TMDICEparticle*> &cur,vector<TMDICEparticle*> &fin);
virtual void stopparticle(TMDICEparticle *p)
{
    stopparticle(p, this->genp,this->genfinp);
};


void make_cascs()
{
    genfin.clear();
    casc.clear();
    vector<TMDICEparticle*>genfintmp, ois;
    genfintmp.clear();
    genfintmp=genfinp;
    ois=root.find(ptrue);
    for(int i=0; i<genfintmp.size();i++)
    {
        genfin.push_back(*genfintmp.at(i));
    }
};

vector<TMDICEparticle*> getcasc()
{
    return root.find(ptrue);
}

vector<TMDICEparticle*> getgenfin()
{
   // return root.find(fincond);
   return genfinp;
}

vector<TMDICEparticle*> gettemp()
{
    return root.find(tempcond);
}

virtual void fillgen(vector<TMDICEparticle*> &cur,vector<TMDICEparticle*> &fin);
virtual void fillgen()
{
    fillgen(genp,genfinp);
};

virtual void initgen(TMDICEparticle &c,vector<TMDICEparticle*> &cur,vector<TMDICEparticle*> &fin);
virtual void initgen()
{
    initgen(root,genp,genfinp);
}

virtual void make_event(TMDICEparticle &c,vector<TMDICEparticle*> &cur,vector<TMDICEparticle*> &fin);

virtual void make_event()
{
    make_event(root,genp,genfinp);
}

};

void determineadress(TMDICEparticle& p1, TMDICEparticle& p2);

#endif
