#include"deps.h"
#include"TMDICE_lib.h"
#include"decor.h"
#include "globalconsts.h"
#include "functions_lib.h"

#include"TMDICE_base.h"
#include"BDIM.h"
#include "medkernels.h"


double probgg;
//map<double,map<double,unordered_map<int,double>>> phiinv;
//map<double,map<double,map<double,double>>> phiinv,difffkt;


map<double,map<double,double>> phmax;
map<double,double> phix;
//const bool simcasc=true;

map<double,map<double,map<double,double>>> phmax_t;
map<double,unordered_map<double,double> > phix_t;
map<double,double> probgg_t;

double cqq[M], cqg[M],cgg[M];
double cqq1[M], cqg1[M],cgg1[M];
double cqq2[M], cqg2[M],cgg2[M];
double cqq3[M], cqg3[M],cgg3[M];
double b1gg,b2gg,bgg;
double b1qg,b2qg,bqg;
double b1qq,b2qq,bqq;
double ccc[3][3][M];

double c_t[3][3][tanzmax][M],c1_t[3][3][tanzmax][M],c2_t[3][3][tanzmax][M],c3_t[3][3][tanzmax][M];
map<double,map<double,map<double,double>>>  b_t,b1_t,b2_t;
//double cqq_u[M], cqg_u[M],cgg_u[M];

//const int tanzmax=1*pow(10,2);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



fkt WW;

map<double,map<double,unordered_map<int,double> >>invph,invph_u;
int anzphiinv=5.*pow(10,3);
int anzphiinvlog=3.*pow(10,2);


void setBDIM()
{
	disclaimer();

	cout<<endl<<"Calculating partition functions: This may take a while."<<endl;
	gevfm=1.0/0.19732;	
	xeps=1.*pow(10,-4);
	xmin=1.*pow(10,-4);
	anz=1.*pow(10,7);

	qmin=.1;
	qmax=40.;	
	qanz=2*pow(10,4);	
	u_min=0.;	
	u_max=pi;
	u_anz=1*pow(10.,7);

	kt0=0.;
	x0=1.;
	typ0=2;
	N_F=3.;
	scatflag=0;
	ktsplitflag=1;
	time_mode=0;
	rotate=1;

	map<string,double*> invars2;
	invars2={{"xeps",&xeps},{"xmin",&xmin},{"qmin",&qmin},{"kt1",&kt0},{"x1",&x0},{"typ1",&typ0},{"nc",&nc},
	{"alphabar",&alphabar},{"alphas",&alphas},{"ndens",&ndens},{"T",&Temp},{"qhat",&qhat},{"emax",&emax},
	{"md",&md},{"taumin",&taumin},{"taumax",&taumax},{"tmin",&tmin},{"tmax",&tmax},{"scat",&scatflag1},{"ktsplit",&ktsplitflag1}
	,{"timemode",&time_mode1},{"rotate",&rotate},{"p1x",&p1x},{"p1y",&p1y},{"p1z",&p1z}};

	setparams(invars2);

	cout<<"time_mode="<<time_mode<<endl;
	cout<<"xmin="<<xmin<<endl;

	ca=nc;
	cf=(nc*nc-1.)/(2.*nc);

	colorfact=pow(cf/ca,1.0);

	if(invars.find("alphabar")==invars.end())
	{
		alphabar=invars["alphas"]/pi;
		alphas=invars["alphas"];
	}
	else
	{
		alphabar=invars["alphabar"];
		alphas=alphabar*pi;
	}
	qhatgev=qhat/gevfm;//inGeV^3
	tstar=alphabar*sqrt(qhatgev/emax);//inGeV
	tstargev=(1.0/tstar);
	tstar=tstargev/gevfm;//infm

	if(invars.find("taumin")!=invars.end() and invars.find("taumax")!=invars.end() ) 
	{
		tmin=taumin*tstar;
		tmax=taumax*tstar;
	}
	if(invars.find("tmin")!=invars.end() and invars.find("tmax")!=invars.end() ) 
	{
		taumin=tmin/tstar;
		taumax=tmax/tstar;
	}

	cout<<"tau="<<tmax/tstar<<" "<<nc*sqrt(nc)*tmax/tstar<<endl;


	if(time_mode==1){Kt=Kt_stat;}
	if(time_mode==2){Kt=Kt_exp;}

	if(scatflag==0){W_scat=0.;WW=WW0;}
	if(scatflag==1){W_scat=4.*pi*nc*alphas*alphas*ndens/(qmin*qmin);W_scat=W_scat*gevfm;WW=WW1;}
	if(scatflag==2){W_scat=(4.*pi*nc*alphas*alphas*ndens/(md*md))*log(1.0+md*md/(qmin*qmin));W_scat=W_scat*gevfm;WW=WW2;}
	if(scatflag==3){md=sqrt(qhatgev/(alphas*nc*Temp)); W_scat=nc*alphas*Temp*pi*qhatgev*log(1.0+md*md/(qmin*qmin));W_scat=W_scat*gevfm;WW=WW2;}

	map<double,map<double,double>> ymap={{1,{{1,5*pow(10,5)}}},{2,{{1,10},{2,1*pow(10,9)}}}};
	map<double,map<double,map<double,double>>> phiinv,difffkt;
	fillgauss(gaussx,gaussw);
	for(double i=2.;i>=1.;i--)
	{
		for(double j=i;j>=1;j--)
		{
			map<double,map<double,string>>pnam={{1,{{1,"qq"},{2,"qg"}}},{2,{{1,"gq"},{2,"gg"}}}};
			cout<<"Calculating "<<pnam[i][j]<<" distributions:"<<endl<<endl;

			fill_partitionfunction_v2(K[i][j],1*xeps,1.0-1*xeps,anz,phiinv[i][j],phmax[i][j]);

			map<double, double>::iterator it,it2;
			map<double,double> m1=phiinv[i][j];
			double m1max=phmax[i][j];

			int anztmp=1;
			invph[i][j][0]=0.;
       		for (it = m1.begin(); it != m1.end(); it++)
			{
				it2=it;
				it2++;
				double ph1tmp=it->first,ph2tmp=it2->first;
				double z1tmp=it->second,z2tmp=it2->second;
				if(anzphiinv*ph1tmp/m1max<anztmp and anzphiinv*ph2tmp/m1max>anztmp )
				{
					invph[i][j][anztmp]=z1tmp+(static_cast<double>(anztmp)*m1max/static_cast<double>(anzphiinv)-ph1tmp)* (z2tmp-z1tmp)/(ph2tmp-ph1tmp);
					anztmp++;
				}
				
			}
			phiinv[i][j].clear();

			double dt=(tmax-tmin)/tanz;

			loadstart=false;
			if(time_mode!=0)
			{
				cout<<"Hi"<<endl;
				for(double t=tmin;t<=tmax;t=t+dt)
				{
					double xx;
					double tnew=tmin+ceil((t-tmin)/dt)*dt;

					map<double,map<double,map<double,map<double,double>>>> phiinv_t,difffkt_t;

					fill_partitionfunction_v2(Kt[i][j],1*xeps,1.0-1*xeps,anz,phiinv_t[i][j][tnew],phmax_t[i][j][tnew],xx,tnew);


					cheb_coeff_split(Kt,phmax_t,phiinv_t,i,j, ymap[i][j],c1_t,c2_t,c3_t,b1_t,b2_t,b_t,tnew);


					if(phmax_t[1.][1.].find(tnew)!=phmax_t[1.][1.].end())
					{
						phix_t[1.][tnew]=phmax_t[1.][1.][tnew];					
					}
					if(phmax_t[2.][2.].find(tnew)!=phmax_t[2.][2.].end() and phmax_t[2.][1.].find(tnew)!=phmax_t[2.][1.].end())
					{
						phix_t[2.][tnew]=phmax_t[2.][2.][tnew]+phmax_t[2.][1.][tnew];
						probgg_t[tnew]=phmax_t[2.][2.][tnew]/phix_t[2.][tnew];
					}
					loadingbar(tmin,tmax,t,tanz);
				}
			}
			
				
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

				cout<<"min="<<exp(m1min_u)<<"max="<<m1max_u<<"delta"<<logdu<<endl;

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
						invph_u[i][j][anztmp_u]=z1tmp+(static_cast<double>(anztmp_u)*m1max_u/static_cast<double>(anzphiinvlog)-ph1tmp)* (z2tmp-z1tmp)/(ph2tmp-ph1tmp);
						cout<<anztmp_u<<" "<<ph1tmp<<" "<<ph2tmp<<" "<<z1tmp<<endl;
						cout<<exp(log(m1max_u)+m1min_u+anztmp_u*logdu)<<endl;
						anztmp_u++;
						cout<<exp(log(m1max_u)+m1min_u+anztmp_u*logdu)<<endl;
					}
					
				}
				phi_u_inv[i][j].clear();
			}
			cout<<endl;
		}
	}
	phix[1.]=phmax[1.][1.];
	phix[2.]=phmax[2.][2.]+phmax[2.][1.];
	probgg=phmax[2.][2.]/phix[2.];


	srand (time(NULL));
	cout<<"Initialization finished"<<endl;
}

double interpolate_phi(double z, double i, double j)
{
	double zind=(z/phmax[i][j])*anzphiinv;
	double zmintmp=max(0.,floor(zind)),zmaxtmp=min(static_cast<double>(anzphiinv),ceil(zind));

	double phmintmp=max(0.,invph[i][j][static_cast<int>(zmintmp)]),phmaxtmp=min(1.,invph[i][j][static_cast<int>(zmaxtmp)]);


	return phmintmp+(zind-zmintmp)*(phmaxtmp-phmintmp);
}

double interpolate_phi_u(double z, double i, double j)
{
				double m1min_u=log(1./anzphiinvlog);//log((ittest->first)/(phi_u_max[i][j]));
				double logdu=-m1min_u/anzphiinvlog;
	//										double utestlog=log(m1max_u)+m1min_u+anztmp_u*logdu;

	double zind=(log(z/phi_u_max[i][j])-m1min_u)/logdu;
	//double zind=((z)/(phi_u_max[i][j]))*anzphiinv;
	//double zmintmp=floor(zind),zmaxtmp=min(static_cast<double>(anzphiinv),ceil(zind));
	double zmintmp=max(0.,floor(zind)),zmaxtmp=min(static_cast<double>(anzphiinvlog),ceil(zind));
	double phmintmp=invph_u[i][j][static_cast<int>(zmintmp)],phmaxtmp=invph_u[i][j][static_cast<int>(zmaxtmp)];

//	double phmintmp=max(0.,invph_u[i][j][static_cast<int>(zmintmp)]),phmaxtmp=min(static_cast<double>(pi),invph_u[i][j][static_cast<int>(zmaxtmp)]);

	//cout<<zmintmp<<" "<<zind<<" "<<zmaxtmp<<" "<<phmintmp<<" "<<phmaxtmp<<endl;
	//cout<<phmintmp+(zind-zmintmp)*(phmaxtmp-phmintmp)<<endl;

	return phmintmp+(zind-zmintmp)*(phmaxtmp-phmintmp);//(zmaxtmp-zmintmp);
	//return invph_u[i][j][static_cast<int>(zind)];
}

double phqq(double z){if(phmax[1][1]>0){return interpolate_phi(z, 1,1);}else{return 0.;}}
double phqg(double z){if(phmax[2][1]>0){return interpolate_phi(z, 2,1);}else{return 0.;}}
double phgg(double z){if(phmax[2][2]>0){return interpolate_phi(z, 2,2);}else{return 0.;}}
map<double,map<double,fkt>>phiinv1={{1.0, {{1.0,phqq}}},{-1.0, {{-1.0,phqq}}},{2.0, {{1.0,phqg},{2.0,phgg}}}};

double phqq(double z,double t){double dt=(tmax-tmin)/tanz;int tind=static_cast<int>((t-tmin)/dt);if(phmax_t[1][1].lower_bound(t)->second>0){return approx_f(z,c1_t[1][1][tind],c2_t[1][1][tind],c3_t[1][1][tind],0.,b1_t[1][1].lower_bound(t)->second,b2_t[1][1].lower_bound(t)->second,b_t[1][1].lower_bound(t)->second,M);;}else{return 0.;}}
double phqg(double z,double t){double dt=(tmax-tmin)/tanz;int tind=static_cast<int>((t-tmin)/dt);if(phmax_t[2][1].lower_bound(t)->second>0){return approx_f(z,c1_t[2][1][tind],c2_t[2][1][tind],c3_t[2][1][tind],0.,b1_t[2][1].lower_bound(t)->second,b2_t[2][1].lower_bound(t)->second,b_t[2][1].lower_bound(t)->second,M);;}else{return 0.;}}
double phgg(double z,double t){double dt=(tmax-tmin)/tanz;int tind=static_cast<int>((t-tmin)/dt);if(phmax_t[2][2].lower_bound(t)->second>0){return approx_f(z,c1_t[2][2][tind],c2_t[2][2][tind],c3_t[2][2][tind],0.,b1_t[2][2].lower_bound(t)->second,b2_t[2][2].lower_bound(t)->second,b_t[2][2].lower_bound(t)->second,M);;}else{return 0.;}}
map<double,map<double,fkt2>>phiinv1_t={{1.0, {{1.0,phqq}}},{-1.0, {{-1.0,phqq}}},{2.0, {{1.0,phqg},{2.0,phgg}}}};

double phi_u_invfkt(double z,double i, double j, double (&c)[M])
{
	map<double,double>::iterator it1=phi_u_inv[i][j].begin(),it2=phi_u_inv[i][j].end();
	it2--;
	return approx_f(z, (it1->first),(it2->first),c,ee_u);
}

//double phqq_u(double z){return phi_u_invfkt(z,1.0,1.0,cqq_u);}
double phqq_u(double z){return interpolate_phi_u(z, 1,1);}
double phqg_u(double z){return interpolate_phi_u(z, 2,1);}
double phgg_u(double z){return interpolate_phi_u(z, 2,2);}
map<double,map<double,fkt>>phi_u_inv1={{1.0, {{1.0,phqq_u}}},{-1.0, {{-1.0,phqq_u}}},{2.0, {{1.0,phqg_u},{2.0,phgg_u}}}};

double phixg(double x){return phix[2.0]/sqrt(x);}
double phixq(double x){return phix[1.0]/sqrt(x);}
map<double,fkt> phi={{1.0,phixq},{2.0,phixg}};
fvec phi2={phixq,phixg};

double approx_t(double tin)
{	
	double dt = (tmax-tmin)/tanz;
	double it = ceil((tin-tmin)/dt);
	double tout = tmin+it*dt;
	tout=tmin+ceil((tout-tmin)/dt)*dt;
	return tout;
}

double phixg_t(double x,double t){	return phix_t[2.0][approx_t(t)]/sqrt(x);}
double phixq_t(double x,double t){ return phix_t[1.0][approx_t(t)]/sqrt(x);}
map<double,fkt2> phi_t={{1.0,phixq_t},{2.0,phixg_t}};

BDIMevent::BDIMevent()
{
	setevolmin(tmax);
	setsetpval(set_to_t);
	setstopcond(stop_with_time);
	setdumpcond(ltxmin);
}

double BDIMevent::evol_select(TMDICEparticle p)
{
	double f=p.typ,x=p.x,told=p.t_old;
	double colorfact2=1.;
	if(abs(f)==1){colorfact2=colorfact;}

	double t1=told,t2;
	double r1,r2,wwww;
	fkt *phitmp=&phi[abs(f)];
	double phi1=(*phitmp)(x);
	if(time_mode!=0)
	{
		do
		{
			r1=rand01();
			t2=(t1)-log(r1)/(up_lim_phi_time[time_mode]*phi[abs(f)](x)/(tstar)+colorfact2*W_scat);

			double dt = (tmax-tmin)/tanz;
			fkt2 *phimtmp=&phi_t[abs(f)];
			double ttmp=max(t2,tmax);
			double phim2=(*phimtmp)(x,ttmp);
			double phim1=(*phimtmp)(x,ttmp-dt);
			double phim=phim1+(ttmp-approx_t(ttmp-dt))*(phim2-phim1)/dt;		
			wwww=(phim/(tstar)+colorfact2*W_scat)/(up_lim_phi_time[time_mode]*phi1/(tstar)+colorfact2*W_scat);

			t1=t2;
			r2=rand01();
		}while(wwww<r2 and t2<tmax);
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
		fkt2 *phitmp=&phi_t[abs(f)];
		phi1=(*phitmp)(x,approx_t(t));
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
		if(r<probgg){f2=2.0;}else{double r2=rand01(); if(r2<0.5){f2=1.0;}else{f2=-1.0;} }
	}
	return f2;
}

double BDIMevent::xselect(TMDICEparticle p,TMDICEparticle p1,TMDICEparticle p2)
{
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
		map<double,double>::iterator it=phmax_t[abs(f1)][abs(f2)].lower_bound(t);
		it--;
		xtmp=select(phiinv1_t[abs(f1)][abs(f2)],it->second ,1*xeps,1.-1*xeps,t);
	}


	return xtmp;
}

double f_u(double u)
{
	return 0.5*(1.-exp(-u)*(sin(u)+cos(u)));
}
void BDIMevent::vertex(TMDICEparticle p,TMDICEparticle &p1,TMDICEparticle &p2)
{
	p1.typ=typselect(p);
	
	p2.typ=typ2[p.typ][p1.typ];

	p1.x=xselect( p, p1, p2);
	
	p2.x=1.0-p1.x;
	p1.kt=p.kt;
	p2.kt=p1.kt;
	
	p1.kt=p1.kt*p1.x;
	p2.kt=p2.kt*p2.x;

	p1.islead=p.islead;
	p2.islead=p.islead;
	if(p1.x<p2.x){p1.islead=0;}
	if(p1.x>p2.x){p2.islead=0;}

	if(ktsplitflag==1)
	{
		double u_tmp1=0,u_tmp2=0;
		
		if(phi_u_max[abs(p.typ)][abs(p1.typ)]>0){u_tmp1=select(phi_u_inv1[abs(p.typ)][abs(p1.typ)],phi_u_max[abs(p.typ)][abs(p1.typ)] ,u_min,u_max);}
	

		double Q_1=2.0*sqrt(qhatgev*p.x*emax/nc)*sqrt(p1.x*(1.0-p1.x)*f[abs(p.typ)][abs(p1.typ)](p1.x) )*u_tmp1;
		Q_1=sqrt(Q_1); 
		double phi_q1=2.0*pi*rand01();
		double Q_2=Q_1;
		double phi_q2=phi_q2+pi;
		if(phi_q2>1.*pi){phi_q2=phi_q2-2*pi;}
		if(phi_q2<-1.*pi){phi_q2=phi_q2+2*pi;}
		
		p1.phik=p.phik;
		p2.phik=p.phik;
		if(rotate==0)
		{
			double p1_kx=p1.kt*cos(p1.phik)+Q_1*cos(phi_q1),p1_ky=p1.kt*sin(p1.phik)+Q_1*sin(phi_q1);
			double p2_kx=p2.kt*cos(p2.phik)+Q_2*cos(phi_q2),p2_ky=p2.kt*sin(p2.phik)+Q_2*sin(phi_q2);


			p1.p0=p.p0*p1.px;
			p1.px=p1_kx;
			p1.py=p1_ky;
			p1.pz=p.pz*p1.x;

			p2.p0=p.p0*p2.px;
			p2.px=p2_kx;
			p2.py=p2_ky;
			p2.pz=p.pz*p2.x;

			p1.kt=sqrt(p1_kx*p1_kx+p1_ky*p1_ky);
			p1.phik=atan2(p1_ky,p1_kx);

			p2.kt=sqrt(p2_kx*p2_kx+p2_ky*p2_ky);
			p2.phik=atan2(p2_ky,p2_kx);
		}
		else
		{
			p1.kt=Q_1;
			p2.kt=Q_2;

			p1.phik=phi_q1;
			p2.phik=phi_q2;

			double ax1_x,ax1_y,ax1_z;
			double ax2_x,ax2_y,ax2_z;
			if(p.pz==0)
			{
				ax1_z=1.;
				ax1_y=0.;
				ax1_x=0.;
			}
			else
			{
				double norm=sqrt(p.py*p.py+p.pz*p.pz);
				ax1_x=0.;
				ax1_y=p1.pz/norm;
				ax1_z=-p1.py/norm;
			}
			ax2_x=ax1_y*p.pz-ax1_z*p.py;
			ax2_y=-(ax1_x*p.pz-ax1_z*p.px);
			ax2_z=ax1_x*p.py-ax2_y*p.px;
			double norm= sqrt(ax2_x*ax2_x+ax2_y*ax2_y+ax2_z*ax2_z);
			ax2_x=ax2_x/norm;
			ax2_y=ax2_y/norm;
			ax2_z=ax2_z/norm;

			p1.p0=p1.x*p.p0;
			p1.px=p1.x*p.px+Q_1*cos(phi_q1)*ax1_x+Q_1*sin(phi_q1)*ax2_x;
			p1.py=p1.x*p.py+Q_1*cos(phi_q1)*ax1_y+Q_1*sin(phi_q1)*ax2_y;
			p1.pz=p1.x*p.pz+Q_1*cos(phi_q1)*ax1_z+Q_1*sin(phi_q1)*ax2_z;

			p2.p0=p2.x*p.p0;
			p2.px=p2.x*p.px+Q_2*cos(phi_q2)*ax1_x+Q_2*sin(phi_q2)*ax2_x;
			p2.py=p2.x*p.py+Q_2*cos(phi_q1)*ax1_y+Q_2*sin(phi_q2)*ax2_y;
			p2.pz=p2.x*p.pz+Q_2*cos(phi_q1)*ax1_z+Q_2*sin(phi_q2)*ax2_z;

			p1.phik=atan2(p1.py,p1.px);
			p1.kt=sqrt(p1.px*p1.px+p1.py*p1.py);
			p2.phik=atan2(p2.py,p2.px);
			p2.kt=sqrt(p2.px*p2.px+p2.py*p2.py);

		}
		
	}
	else
	{
		p1.phik=p.phik;
		p2.phik=p.phik;

		p1.p0=p1.x*p.p0;
		p1.px=p1.x*p.px;
		p1.py=p1.x*p.py;
		p1.pz=p1.x*p.pz;

		p2.p0=p2.x*p.p0;
		p2.px=p2.x*p.px;
		p2.py=p2.x*p.py;
		p2.pz=p2.x*p.pz;		
	}
	
	p1.x=p1.x*p.x;
	p2.x=p2.x*p.x;
	
	p1.t_old=p.t;
	p2.t_old=p.t;
}

void BDIMevent::scat(TMDICEparticle p, TMDICEparticle &p1)
{
	double ktnew=WW(rand01());
	ktnew=sqrt(ktnew);
	double phinew=2.*pi*rand01();

		p1.islead=p.islead;


	double ktold=p1.kt;
	double p1x=p.kt*cos(p.phik)+ktnew*cos(phinew),p1y=p.kt*sin(p.phik)+ktnew*sin(phinew);

	p1.phik=atan2(p1y,p1x);
	p1.kt=sqrt(p1x*p1x+p1y*p1y);

	if(rotate==0)
	{
		p1.px=p1x;
		p1.py=p1y;
	}
	else
	{
		double ax1_x,ax1_y,ax1_z;
		double ax2_x,ax2_y,ax2_z;
		if(p.pz==0)
		{
			ax1_z=1.;
			ax1_y=0.;
			ax1_x=0.;
		}
		else
		{
			double norm=sqrt(p.py*p.py+p.pz*p.pz);
			ax1_x=0.;
			ax1_y=p1.pz/norm;
			ax1_z=-p1.py/norm;
		}
		ax2_x=ax1_y*p.pz-ax1_z*p.py;
		ax2_y=-(ax1_x*p.pz-ax1_z*p.px);
		ax2_z=ax1_x*p.py-ax2_y*p.px;
		double norm= sqrt(ax2_x*ax2_x+ax2_y*ax2_y+ax2_z*ax2_z);
		ax2_x=ax2_x/norm;
		ax2_y=ax2_y/norm;
		ax2_z=ax2_z/norm;

		p1.p0=p.p0;
		p1.px=p.px+p1.kt*cos(p1.phik)*ax1_x+p1.kt*sin(p1.phik)*ax2_x;
		p1.py=p.py+p1.kt*cos(p1.phik)*ax1_y+p1.kt*sin(p1.phik)*ax2_y;
		p1.pz=p.pz+p1.kt*cos(p1.phik)*ax1_z+p1.kt*sin(p1.phik)*ax2_z;
		
		p1.phik=atan2(p1.py,p1.px);
		p1.kt=sqrt(p1.px*p1.px+p1.py*p1.py);
	}
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

void  BDIMevent::fillgen()
{
	vector<TMDICEparticle> newgen,newgenfin;
	if(gen.size()>0)
	{
		for( int j=0; j<gen.size(); j++)
		{
			TMDICEparticle p=gen.at(j);
			
			string sc=interact_select(p);
			
			if(sc=="b")
			{
				TMDICEparticle p1,p2;
				vertex(p,p1,p2);
				p1.t=evol_select(p1);
				p2.t=evol_select(p2);
				
				stopparticle(p1,newgen,newgenfin);
				stopparticle(p2,newgen,newgenfin);
			}
			
			if(sc=="c")
			{
				TMDICEparticle p1=p;
				scat(p,p1);

				p1.t_old=p1.t;
				p1.t=evol_select(p1);
				
				stopparticle(p1,newgen,newgenfin);
			}
		}
	
	gen=newgen;
	genfin.insert(genfin.end(),newgenfin.begin(),newgenfin.end());
	casc.insert(casc.end(),gen.begin(),gen.end());
	casc.insert(casc.end(),newgenfin.begin(),newgenfin.end());
	}
}