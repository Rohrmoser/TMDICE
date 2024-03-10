#include"deps.h"
#include"TMDICE_lib.h"
#include"decor.h"
#include "globalconsts.h"
#include "functions_lib.h"

#include"TMDICE_base.h"
#include"BDIM.h"
#include "medkernels.h"

using namespace std;

double probgg;

map<double,map<double,double>> phmax;
map<double,double> phix;

map<double,map<double,map<double,double>>> phmax_t;
map<double,map<double,double> > phix_t,Delta_t,iDelta_t;
map<double,double> probgg_t;

fkt WW;

map<double,map<double,unordered_map<int,double> >>/*invph,*/invph_u;
map<double,map<double,map<double,double> >>invph;

int anzphiinvlog=3.*pow(10,2);
map<double,map<double,map<double,map<double,double>>>> invph_t;

template<typename ...Rest> void fillphi(double K(double xx, Rest...rest),std::map<double,double>&invphi,double &phimax,double x1, 
double x2/*,double xx,*/,double dx,Rest...rest)
{
	double w=0;

	double x=x1;
	//double dx=0.01;
	while(x<x2)
	{
		double dxx=std::min(dx,1.-xeps-x),fac=0.1;//fac=0.02;
		while(abs(K(x,rest...)-K(x+dxx,rest...))>fac*K(x,rest...)){dxx*=fac;}
		w=w+(K(x,rest...)+K(x+dxx,rest...))*0.5*dxx;
		invphi[w]=x+dxx;
		x+=dxx;         
	}
	double dxx=1.-xeps-x;//fac=0.02;
	w=w+(K(x,rest...)+K(x+dxx,rest...))*0.5*dxx;
	invphi[w]=x+dxx;

	phimax=w;
	//cout<<floor(100.*t/tmax)<<"% calculated"<<'\r';
}

void fillphi_t(double i, double j)
{
	double dt=0.001;//(tmax-tmin)/tanz;
	//double dx=0.01;
    for(double t=tmin;t<=2*tmax;t+=dt)
    {
		double xx;
		fillphi(Kt[i][j],invph_t[i][j][t],phmax_t[i][j][t],xeps, 1.-xeps, /*xx,*/0.01,t);
		cout<<floor(100.*t/(2*tmax))<<"% calculated"<<'\r';
    }
    cout<<"calculated       \n";
}

void setkernelsbdim()
{
	if(time_mode==1){Kt=Kt_stat;}
	if(time_mode==2){Kt=Kt_exp;}
	if(scatflag==0){W_scat=0.;WW=WW0;}
	if(scatflag==1){W_scat=4.*pi*nc*alphas*alphas*ndens/(qmin*qmin);W_scat=W_scat*gevfm;WW=WW1;}
	if(scatflag==2){W_scat=(4.*pi*nc*alphas*alphas*ndens/(md*md))*log(1.0+md*md/(qmin*qmin));W_scat=W_scat*gevfm;WW=WW2;}
	if(scatflag==3){md=sqrt(qhatgev/(alphas*nc*Temp)); W_scat=nc*alphas*Temp*pi*qhatgev*log(1.0+md*md/(qmin*qmin));W_scat=W_scat*gevfm;WW=WW2;}
}

void obtaincdfsbdim()
{
	//map<double,map<double,map<double,double>>> phiinv;
	for(double i=2.;i>=1.;i--)
	{
		for(double j=i;j>=1;j--)
		{
			map<double,map<double,string>>pnam={{1,{{1,"qq"},{2,"qg"}}},{2,{{1,"gq"},{2,"gg"}}}};
			cout<<"Calculating "<<pnam[i][j]<<" distributions:"<<endl<<endl;

			if(time_mode!=0){fillphi_t(i,j);}


			if(time_mode==0)
			{fillphi(K[i][j],invph[i][j],phmax[i][j],xeps,1.0-xeps,0.001);}
			//double dt=(tmax-tmin)/tanz;
		
			if(ktsplitflag==1)
			{
				fill_partitionfunction_v2(K_u[i][j],u_min,pi,u_anz,phi_u_inv[i][j],phi_u_max[i][j]);

				map<double, double>::iterator it_u,it2_u;
				map<double,double> m1_u=phi_u_inv[i][j];
				double m1max_u=(phi_u_max[i][j]);
				map<double,double>::iterator ittest=m1_u.begin();
				ittest++;
				double m1min_u=log(1./anzphiinvlog);
				double logdu=-m1min_u/anzphiinvlog;

				int anztmp_u=1;
				invph_u[i][j][0]=0.;
				for (it_u = m1_u.begin(); it_u != m1_u.end(); it_u++)
				{
					double utestlog=log(m1max_u)+m1min_u+anztmp_u*logdu;
					
					it2_u=it_u;
					it2_u++;
					double ph1tmp=(it_u->first),ph2tmp=(it2_u->first);
					double z1tmp=it_u->second,z2tmp=it2_u->second;
					if(exp(utestlog)>=ph1tmp and exp(utestlog)<ph2tmp)
					{
						invph_u[i][j][anztmp_u]=z1tmp+(static_cast<double>(anztmp_u)*m1max_u
						/static_cast<double>(anzphiinvlog)-ph1tmp)* (z2tmp-z1tmp)/(ph2tmp-ph1tmp);
						anztmp_u++;
					}
					
				}
				phi_u_inv[i][j].clear();
			}
		}
	}
	phix[1.]=phmax[1.][1.];
	phix[2.]=phmax[2.][2.]+phmax[2.][1.];
	probgg=phmax[2.][2.]/phix[2.];


	cout<<"starting time:"<<tmin<<"\n";
	Delta_t[1.][tmin]=0.;
	Delta_t[2.][tmin]=0.;
	iDelta_t[1.][0.]=tmin;
	iDelta_t[2.][0.]=tmin;

	double dt=0.001,wqold=0.,wgold=0.;

	map<double,double>::iterator it;
	for(it=phmax_t[1.][1.].begin();it!=phmax_t[1.][1.].end();it++)
	{
		double t=it->first;
		phix_t[1.][t]=(phmax_t[1.][1.].lower_bound(t))->second;
		phix_t[2.][t]=(phmax_t[2.][1.].lower_bound(t))->second+(phmax_t[2.][2.].lower_bound(t))->second;
		probgg_t[t]=(phmax_t[2.][2.].lower_bound(t))->second/(phix_t[2.].lower_bound(t))->second;
		if(t!=tmin)
		{
			Delta_t[1.][t]=wqold+(0.5*dt*((phix_t[1.].lower_bound(t))->second+(phix_t[1.].lower_bound(t-dt))->second));
			wqold=(Delta_t[1.].lower_bound(t))->second;
			iDelta_t[1.][wqold]=t;

			Delta_t[2.][t]=wgold+(0.5*dt*((phix_t[2.].lower_bound(t))->second+(phix_t[2.].lower_bound(t-dt))->second));
			wgold=(Delta_t[2.].lower_bound(t))->second;
			iDelta_t[2.][wgold]=t;
		}
	}
}

double interpolate(double w, map<double,double>grid)
{
    map<double,double>::iterator it1=grid.lower_bound(w), 
    it2=grid.upper_bound(w);
    it1--;
    double wm=it1->first;
    double wu=it2->first;
    double xm=it1->second;
    double xu=it2->second;
    return xm+(xu-xm)*(w-wm)/(wu-wm);
}

double interpolate_phi_u(double z, double i, double j)
{
	double m1min_u=log(1./anzphiinvlog);
	double logdu=-m1min_u/anzphiinvlog;
	double zind=(log(z/phi_u_max[i][j])-m1min_u)/logdu;
	double zmintmp=max(0.,floor(zind)),zmaxtmp=min(static_cast<double>(anzphiinvlog),ceil(zind));
	double phmintmp=invph_u[i][j][static_cast<int>(zmintmp)],phmaxtmp=invph_u[i][j][static_cast<int>(zmaxtmp)];
	return phmintmp+(zind-zmintmp)*(phmaxtmp-phmintmp);
}

double phqq(double z){if(phmax[1][1]>0){return interpolate(z,invph[1.][1.]);}else{return 0.;}}
double phqg(double z){if(phmax[2][1]>0){return interpolate(z,invph[2.][1.]);}else{return 0.;}}
double phgg(double z){if(phmax[2][2]>0){return interpolate(z,invph[2.][2.]);}else{return 0.;}}
map<double,map<double,fkt>>phiinv1={{1.0, {{1.0,phqq}}},{-1.0, {{-1.0,phqq}}},{2.0, {{1.0,phqg},{2.0,phgg}}}};

double select(fkt2 Finv, double Fmax/*,double wmin,double wmax*/,double t)
{
	double r=Fmax*rand01();
	double w=Finv(r,t);
	return w;
}


double phqq(double z, double t){if((phmax_t[1][1].lower_bound(t)--)->second>0){return interpolate(z,(invph_t[1.][1.].lower_bound(t)--)->second);}else{return 0.;}}
double phqg(double z, double t){if((phmax_t[2][1].lower_bound(t)--)->second>0){return interpolate(z,(invph_t[2.][1.].lower_bound(t)--)->second);}else{return 0.;}}
double phgg(double z,double t){if((phmax_t[2][2].lower_bound(t)--)->second>0){return interpolate(z,(invph_t[2.][2.].lower_bound(t)--)->second);}else{return 0.;}}
map<double,map<double,fkt2>>phiinv1_t={{1.0, {{1.0,phqq}}},{-1.0, {{-1.0,phqq}}},{2.0, {{1.0,phqg},{2.0,phgg}}}};

double phqq_u(double z){return interpolate_phi_u(z, 1,1);}
double phqg_u(double z){return interpolate_phi_u(z, 2,1);}
double phgg_u(double z){return interpolate_phi_u(z, 2,2);}
map<double,map<double,fkt>>phi_u_inv1={{1.0, {{1.0,phqq_u}}},{-1.0, {{-1.0,phqq_u}}},{2.0, {{1.0,phqg_u},{2.0,phgg_u}}}};

double phixg(double x){return phix[2.0]/sqrt(x);}
double phixq(double x){return phix[1.0]/sqrt(x);}
map<double,fkt> phi={{1.0,phixq},{2.0,phixg}};

double phixg_t(double x,double t){	return (phix_t[2.0].lower_bound(t)--)->second/sqrt(x);}
double phixq_t(double x,double t){ return (phix_t[1.0].lower_bound(t)--)->second/sqrt(x);}
map<double,fkt2> phi_t={{1.0,phixq_t},{2.0,phixg_t}};

BDIMevent::BDIMevent()
{
	setevolmin(tmax);
	setsetpval(set_to_t);
	setstopcond(stop_with_time);
	setdumpcond(ltxmin);
	cascfill=true;
}

double BDIMevent::evol_select(TMDICEparticle p)
{
	double f=p.typ,x=p.x,told=p.old_part->t;
	double colorfact2=1.;
	if(abs(f)==1){colorfact2=colorfact;}

	double t1=told,t2;
	double r1,r2,wwww;
	fkt *phitmp=&phi[abs(f)];
	double phi1=(*phitmp)(x);
	if(time_mode!=0)
	{
		double r1=rand01();
		t2=iDelta_t[abs(f)].lower_bound(-sqrt(x)*tstar*log(r1)+Delta_t[abs(f)].lower_bound(t1)->second)->second;
	}
	if(time_mode==0)
	{
		double r=rand01(),tfin;
		t2=(t1)-log(r)/(phi1/(tstar)+colorfact2*W_scat); 
	}
	return t2;
}

string BDIMevent::interact_select(TMDICEparticle p)
{
	double f=round(abs(p.typ)),x=p.x,t=p.t;
	string itp;
	double phi1;
	if(time_mode==0)
	{
		fkt *phitmp=&phi[abs(f)];
		phi1=(*phitmp)(x);
	}
	if(time_mode!=0)
	{
		//the time dependent case is collinear:
		itp="b";
	}	
	double colorfact2=1.;
	if(abs(f)==1){colorfact2=pow(colorfact,1);}
	double r=rand01();
	if(r<=phi1/(phi1+tstar*colorfact2*W_scat) ){itp="b";}
	if(r>phi1/(phi1+tstar*colorfact2*W_scat) ){itp="c";}
	return itp;
}

double BDIMevent::typselect(TMDICEparticle p)
{
	double f1=round(p.typ),t=p.t,x=p.x;
	double f2;
	double r=rand01();
	if(f1==1.0){f2=1.0;}
	if(f1==-1.0){f2=-1.0;}

	if(f1==2.0)
	{
		double rcp;
		if(time_mode==0)
		{
			rcp=probgg;
		}
		
		if(time_mode!=0)
		{
			map<double,double> *probgg_t1=&probgg_t;
			map<double,double>::iterator it=(probgg_t1->lower_bound(t) );
			(it)--;
			rcp=it->second;
		}
		if(r<rcp){f2=2.0;}else{double r2=rand01(); if(r2<0.5){f2=1.0;}else{f2=-1.0;} }

	}
	return f2;
}

double BDIMevent::xselect(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
	//cout<<"x"<<"\n";
	double f1=round(p.typ);
	double f2=round(p1.typ);
	double t=p.t;
	double xtmp;
	double r=rand01()*anz*dinorm;
	unordered_map<int,double> *phinvtmp;

	if(time_mode==0)
	{
		xtmp=select(phiinv1[abs(f1)][abs(f2)],phmax[abs(f1)][abs(f2)] ,1*xeps,1.-1*xeps);

	}
	if(time_mode!=0)
	{

		t=min(t,tmax);
		xtmp=select(phiinv1_t[abs(f1)][abs(f2)],(phmax_t[abs(f1)][abs(f2)].lower_bound(t)--)->second,t);
	}
	return xtmp;
}

void BDIMevent::vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2)
{
	//cout<<"vertex"<<endl;
	p1=p;
	p2=p;
	p1.typ=typselect(p);
	
	p2.typ=typ2[p.typ][p1.typ];
	p1.x=xselect( p, p1, p2);
	p2.x=1.0-p1.x;

	
	p1.kt=p.kt*p1.x;
	p2.kt=p.kt*p2.x;

	if(p1.x<p2.x){p1.islead=0;p1.adr+="2";p2.adr+="1";}
	else{p2.islead=0;p2.adr+="2";p1.adr+="1";}

	if(ktsplitflag==1)
	{
		double u_tmp1=0,u_tmp2=0;
		
		if(phi_u_max[abs(p.typ)][abs(p1.typ)]>0){u_tmp1=select(phi_u_inv1[abs(p.typ)][abs(p1.typ)],phi_u_max[abs(p.typ)][abs(p1.typ)] ,u_min,u_max);}
	

		double Q_1=2.0*sqrt(qhatgev*p.x*emax/nc)*sqrt(p1.x*(1.0-p1.x)*f[abs(p.typ)][abs(p1.typ)](p1.x) )*u_tmp1;
		Q_1=sqrt(Q_1); 
		double phi_q1=2.0*pi*rand01();
		double Q_2=Q_1;
		double phi_q2=phi_q1+pi;
		if(phi_q2>1.*pi){phi_q2=phi_q2-2*pi;}
		if(phi_q2<-1.*pi){phi_q2=phi_q2+2*pi;}

		if(rotate==0)
		{
			double p1_kx=p1.kt*cos(p1.phik)+Q_1*cos(phi_q1),p1_ky=p1.kt*sin(p1.phik)+Q_1*sin(phi_q1);
			double p2_kx=p2.kt*cos(p2.phik)+Q_2*cos(phi_q2),p2_ky=p2.kt*sin(p2.phik)+Q_2*sin(phi_q2);

			fourmom tmpmom1(sqrt(p1_kx*p1_kx+p1_ky*p1_ky+p1.x*p.p.pz*p1.x*p.p.pz),p1_kx,p1_ky,p1.x*p.p.pz);
			fourmom tmpmom2(sqrt(p2_kx*p2_kx+p2_ky*p2_ky+p2.x*p.p.pz*p2.x*p.p.pz),p2_kx,p2_ky,p2.x*p.p.pz);

			p1.p=tmpmom1;
			p2.p=tmpmom2;

			p1.kt=sqrt(p1_kx*p1_kx+p1_ky*p1_ky);
			p1.phik=atan2(p1_ky,p1_kx);

			p2.kt=sqrt(p2_kx*p2_kx+p2_ky*p2_ky);
			p2.phik=atan2(p2_ky,p2_kx);
		}
		else
		{
			double p1_kx=Q_1*cos(phi_q1),p1_ky=Q_1*sin(phi_q1);
			double p2_kx=Q_2*cos(phi_q2),p2_ky=Q_2*sin(phi_q2);
			fourmom p1oldframe(p1.x*sqrt(p.p.sk3(p.p)),p1_kx,p1_ky,p1.x*sqrt(p.p.sk3(p.p)));
			fourmom p2oldframe(p2.x*sqrt(p.p.sk3(p.p)),p2_kx,p2_ky,p2.x*sqrt(p.p.sk3(p.p)));
			fourmom tmp1=Rotate(p.p,p1oldframe);
			p1.p=tmp1;
			fourmom tmp2=Rotate(p.p,p2oldframe);
			p2.p=tmp2;	

			p1.kt=sqrt(p1_kx*p1_kx+p1_ky*p1_ky);
			p1.phik=atan2(p1_ky,p1_kx);

			p2.kt=sqrt(p2_kx*p2_kx+p2_ky*p2_ky);
			p2.phik=atan2(p2_ky,p2_kx);
		}
		
	}
	else
	{
		fourmom tmpmom1(p1.x*p.p.p0,p1.x*p.p.px,p1.x*p.p.py,p1.x*p.p.pz);
		fourmom tmpmom2(p2.x*p.p.p0,p2.x*p.p.px,p2.x*p.p.py,p2.x*p.p.pz);

		p1.p=tmpmom1;
		p2.p=tmpmom2;
	}
	
	p1.x=p1.x*p.x;
	p2.x=p2.x*p.x;

	p1.old_part=&p;
	p2.old_part=&p;


}

void BDIMevent::scat(TMDICEparticle p, TMDICEparticle &p1)
{
	//cout<<"scat"<<endl;
	double ktnew=WW(rand01());
	ktnew=sqrt(ktnew);
	double phinew=2.*pi*rand01();

	double ktold=p1.kt;
	double px1=p.kt*cos(p.phik)+ktnew*cos(phinew),py1=p.kt*sin(p.phik)+ktnew*sin(phinew);

	p1.phik=atan2(py1,px1);
	p1.kt=sqrt(px1*px1+py1*py1);

	if(rotate==0)
	{
		fourmom tmpmom(sqrt(px1*px1+py1*py1+p.p.pz*p.p.pz),px1,py1,p.p.pz);
		p1.p=tmpmom;
	}
	else
	{
		px1=ktnew*cos(phinew);
		py1=ktnew*sin(phinew);

		p1.phik=atan2(py1,px1);
		p1.kt=sqrt(px1*px1+py1*py1);
		fourmom p1oldframe(p.p.p0,px1,py1,p.p.pz);
		fourmom tmp1=Rotate(p.p,p1oldframe);

		p1.p=tmp1;
	}

	p1.old_part=&p;

}

//dumpcond
bool ltxmin(TMDICEparticle p)
{
	return p.x<=xmin;
}
//stopcond
bool stop_with_time(TMDICEparticle p)
{
	return p.t>=tmax;
}
//setpval
void set_to_t(TMDICEparticle &p, double tt)
{
	p.t=tt;
}

void  BDIMevent::fillgen(vector<TMDICEparticle*> &cur,vector<TMDICEparticle*> &fin)
{
	vector<TMDICEparticle*>genptmp=cur;
	cur.clear();
	for( int j=0; j<genptmp.size(); j++)
	{
		TMDICEparticle* p=genptmp.at(j);
		
		string sc=interact_select(*p);
		
		if(sc=="b")
		{
			TMDICEparticle p1,p2;
			vertex(*p,p1,p2);

			p->insert(p1);
			p->insert(p2);
			p->left->t=evol_select(*(p->left));
			p->right->t=evol_select(*(p->right));

			stopparticle(p->left,cur,fin);
			stopparticle(p->right,cur,fin);
		}			
		if(sc=="c")
		{
			TMDICEparticle p1=*p;
			scat(*p,p1);
			p->insert(p1);

			p->left->t=evol_select(*(p->left));

			stopparticle(p->left,cur,fin);
		}
	}
}