#include"deps.h"

#include"TMDICE_lib.h"
#include"TMDICE.h"

double gevfm;

double xeps,xmin;
double anz;
double qmin;
double qmax;
double qanz;
double u_min;
double u_max;
double u_anz;

double colorfact,ca,cf;


bool func_exec_tmdice=false;

void disclaimer()
{
	if(func_exec_tmdice==false)
	{
	cout<<endl;
	cout<<"This is the TMDICE program"<<endl;
	cout<<"by Martin Rohrmoser, 2021."<<endl;
	cout<<"If you use it, please cite it!"<<endl<<endl;
	}
	func_exec_tmdice=true;
}

double tmin,tmax,taumin,taumax, qhat,qhatgev, emax, alphas, alphabar,nc, ndens, inittyp,W_scat, md,tstar,tstargev,Temp,x0,kt0,typ0; 
int ktsplitflag, scatflag;

map<string,double> invars;


double probgg;
map<double,map<double,map<double,double>>> phiinv;
map<double,map<double,double>> phmax;
map<double,double> phix;
const bool simcasc=true;

double kqq(double z){return 1*0.5*cf*(1.0+z*z)*pow(1.0-z,-1)*sqrt(z*ca+pow(1.0-z,2)*cf)/sqrt(z*(1.0-z));}
double kqg(double z){double wert=1*N_F*T_R*(z*z+pow(1.0-z,2)) *sqrt(cf-z*(1.0-z)*ca)/sqrt(z*(1.0-z));return wert/(2.);}
double kgq(double z){return 1*0.5*cf*(1.0+pow(1.0-z,2))*(1.0/z)*sqrt((1.0-z)*ca+z*z*cf)/sqrt(z*(1.0-z));}
double kgg(double z){double wert=pow(ca,1.5)*pow(1.0-z+z*z,2.5)*pow(z*(1.0-z),-1.5); return wert/2.;}
map<double,map<double,fkt>>K={{1.0, {{1.0,kqq},{2.0,kgq}}},{-1.0, {{-1.0,kqq},{2.0,kgq}}},{2.0, {{1.0,kqg},{2.0,kgg}}}};

double k_u_qq(double z){return sin(z)*exp(-z);;}
double k_u_qg(double z){return sin(z)*exp(-z);;}
double k_u_gq(double z){return sin(z)*exp(-z);;}
double k_u_gg(double z){return sin(z)*exp(-z);}
map<double,map<double,fkt>>K_u={{1.0, {{1.0,k_u_qq},{2.0,k_u_gq}}},{-1.0, {{-1.0,k_u_qq},{2.0,k_u_gq}}},{2.0, {{1.0,k_u_qg},{2.0,k_u_gg}}}};

double fgg(double z){return 1*((1.-z)+z*z)*(ca);}
double fqg(double z){return cf-z*(1.-z)*ca;}
double fgq(double z){return (1.-z)*ca+z*z*cf;}
double fqq(double z){return z*ca+(1.-z)*(1.-z)*cf;}
map<double,map<double,fkt>>f={{1.0, {{1.0,fqq},{2.0,fgq}}},{-1.0, {{-1.0,fqq},{2.0,fgq}}},{2.0, {{1.0,fqg},{2.0,fgg}}}};

map<double,map<double,map<double,double>>> phi_u_inv;
map<double,map<double,double>> phi_u_max;

double WW0(double R){return 0;}
double WW1(double R){return qmin*qmin/(1.0-R);}//+
double WW2(double R){return md*md/(pow(1.0+md*md/(qmin*qmin),1.0-R)-1.0);}//+
double WW3(double R){return md*md/(pow(0.5*(1.0+md*md/(qmin*qmin) ),1.0-R)-1.0);}//+

fkt WW;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const int M=100;
const double m=static_cast<double>(M);

const double ee=1*pow(10,0),ee_u=7*pow(10,-5);

double cqq[M], cqg[M],cgg[M];
double cqq1[M], cqg1[M],cgg1[M];
double cqq2[M], cqg2[M],cgg2[M];
double cqq3[M], cqg3[M],cgg3[M];
double b1gg,b2gg,bgg;
double b1qg,b2qg,bqg;
double b1qq,b2qq,bqq;
double cqq_u[M], cqg_u[M],cgg_u[M];

vector<double> det_boundaries(fkt f, double y, double zmin,double zmax)
{
	vector<double> t;
	double dz=pow(10,-3);
	for(double z=zmin;z<zmax;z=z+dz )
	{
		if((f(z)>1./y and f(z+dz)<=1./y) or (f(z)<1./y and f(z+dz)>=1./y))
		{
			double zz;
			double ival=integral(f,zmin,z,zz);
			t.push_back(ival);
		}
	}
	return t;
}

void make_c(map<double,double> f,double (&c)[M],double a, double b,double eee)
{
	map<double,double>::iterator it2=f.end(),it_tmp;
	it2--;
	a=a+eee;
	b=b+eee;
	
	for(int j=0;j<=M-1;j++)
	{
		c[j]=0;
		for(double k=1.0;k<=m;k++)
		{
			double x=( 0.5*(a+b +cos(pi*(k-0.5)/m)*(b-a)) )-eee;
			it_tmp=f.lower_bound(x);
			x=it_tmp->first;
			c[j]=c[j]+((f[x]))*cos(pi*static_cast<double>(j)*(k-0.5)/m);
		}
		c[j]=c[j]*2.0/m;
	}
}

void make_3c(map<double,double>f,double (&c1)[M],double (&c2)[M],double (&c3)[M],double a,double b1,double b2,double b)
{
	make_c(f,c1,a,b1,0.);
	make_c(f,c2,b1,b2,0.);
	make_c(f,c3,b2,b,0.);
}

double approx_f(double x, double a, double b,double (&c)[M],double eee)
{
	double y=((x+eee)-0.5*((a+eee)+(b+eee)))/(0.5*((b+eee)-(a+eee)));
	double f=0.5*c[0];
	for(int k=1;k<=M-1;k++)
	{
		f=f+c[k]*T(y,k);
	}
	return (f);
}

double approx_f(double x,double (&c1)[M],double (&c2)[M],double (&c3)[M],double a,double b1,double b2,double b)
{
	double aa,bb,f;
	if(a<b)
	{
		if(x<=b1)
		{
			aa=a,bb=b1;
			double y=((x)-0.5*((aa)+(bb)))/(0.5*((bb)-(aa)));

			f=0.5*c1[0];
			for(int k=1;k<=M-1;k++)
			{
				f=f+c1[k]*T(y,k);
			}
		}
		if(x>b1 and x<=b2)
		{
			aa=b1;bb=b2;
			double y=((x)-0.5*((aa)+(bb)))/(0.5*((bb)-(aa)));
			f=0.5*c2[0];
			for(int k=1;k<=M-1;k++)
			{
				f=f+c2[k]*T(y,k);
			}
		}
		if(x>b2)
		{
			aa=b2;bb=b;
			double y=((x)-0.5*((aa)+(bb)))/(0.5*((bb)-(aa)));
			f=0.5*c3[0];
			for(int k=1;k<=M-1;k++)
			{
				f=f+c3[k]*T(y,k);
			}
		}
	}
	else
	{
		f=0;
	}

	return f;
}



void readTMDICEparameters(string snom)
{
	disclaimer();
	string comment;
	ifstream s;
	cout<<endl<<"This is the name of the input file: "<<snom<<endl;
	s.open(snom.c_str() );
	cout<<endl<<"Input parameters:"<<endl;
	while(s.eof()==false )
	{
		string vnam;
		double vval;
		s>>vnam>>vval;
		getline(s,comment);
		if(vnam!="")
			{
				invars[vnam]=vval;
				cout<<vnam<<" "<<invars[vnam]<<comment<<endl;
			}
	}
	s.close();
	string errmess="The following information is necessary: 1.) nc and ndens, qhat and emax 2.)Either taumin and taumax or tmin and tmax 3.) alphas or alphabar";
	auto fin=invars.end();
	if(invars.find("nc")==fin or invars.find("ndens")==fin or 
	invars.find("qhat")==fin or invars.find("emax")==fin ){cout<<errmess<<endl; exit(1);}
	if((invars.find("taumin")==fin or invars.find("taumax")==fin) and (invars.find("tmin")==fin or invars.find("tmax")==fin)){cout<<errmess<<endl;exit(2);}
	if(invars.find("alphas")==fin and invars.find("alphabar")==fin){cout<<errmess<<endl;exit(3);}
}


void readTMDICEparameters(map<string,double>iv1)
{
	disclaimer();
	
	invars=iv1;
	string errmess="The following information is necessary: 1.) nc and ndens, qhat and emax 2.)Either taumin and taumax or tmin and tmax 3.) alphas or alphabar";
	auto fin=invars.end();
	if(invars.find("nc")==fin or invars.find("ndens")==fin or 
	invars.find("qhat")==fin or invars.find("emax")==fin ){cout<<errmess<<endl; exit(1);}
	if((invars.find("taumin")==fin or invars.find("taumax")==fin) and (invars.find("tmin")==fin or invars.find("tmax")==fin)){cout<<errmess<<endl;exit(2);}
	if(invars.find("alphas")==fin and invars.find("alphabar")==fin){cout<<errmess<<endl;exit(3);}
}

void setTMDICE()
{
	disclaimer();


	cout<<endl<<"Calculating partition functions: This may take a while."<<endl;
	gevfm=1.0/0.19732;
	
	xeps=1.*pow(10,-4);
	if(invars.find("xeps")!=invars.end()){xeps=invars["xeps"];}
	
	xmin=1.*pow(10,-4);
	if(invars.find("xmin")!=invars.end()){xmin=invars["xmin"];}

	anz=1*pow(10,4);

	qmin=.1;
	if(invars.find("qmin")!=invars.end()){qmin=invars["qmin"];}
	
	qmax=40.;
	
	qanz=2*pow(10,4);
	
	u_min=0.;
	
	u_max=pi;

	u_anz=1*pow(10.,4);

	kt0=0.;
	if(invars.find("kt1")!=invars.end()){kt0=invars["kt1"];}
	x0=1.;
	if(invars.find("x1")!=invars.end()){x0=invars["x1"];}
	typ0=2;
	if(invars.find("typ1")!=invars.end()){typ0=invars["typ1"];}

	
	N_F=3.;
	nc=invars["nc"];

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

	ndens=invars["ndens"];//in GeV^3
	if(invars.find("T")!=invars.end() )
	{
		Temp=invars["T"];//in GeV
	}
	qhat=invars["qhat"];
	qhatgev=invars["qhat"]/gevfm;//inGeV^3
	emax=invars["emax"];
	md=invars["md"];
	tstar=alphabar*sqrt(qhatgev/emax);//inGeV
	tstargev=(1.0/tstar);
	tstar=tstargev/gevfm;//infm

	if(invars.find("taumin")!=invars.end() and invars.find("taumax")!=invars.end() ) 
	{
		taumin=invars["taumin"];
		taumax=invars["taumax"];
		tmin=taumin*tstar;
		tmax=taumax*tstar;
	}
	if(invars.find("tmin")!=invars.end() and invars.find("tmax")!=invars.end() ) 
	{
		tmin=invars["tmin"];
		tmax=invars["tmax"];
		taumin=tmin/tstar;
		taumax=tmax/tstar;
	}

	fillgauss(gaussx,gaussw);

	scatflag=0;
	ktsplitflag=1;
	
	scatflag=static_cast<int>(invars["scat"]);
	ktsplitflag=static_cast<int>(invars["ktsplit"]);

	if(scatflag==0){W_scat=0.;WW=WW0;}
	if(scatflag==1){W_scat=4.*pi*nc*alphas*alphas*ndens/(qmin*qmin);W_scat=W_scat*gevfm;WW=WW1;}
	if(scatflag==2){W_scat=(4.*pi*nc*alphas*alphas*ndens/(md*md))*log(1.0+md*md/(qmin*qmin));W_scat=W_scat*gevfm;WW=WW2;}
	if(scatflag==3){W_scat=nc*alphas*Temp*log(1.0+md*md/(qmin*qmin));W_scat=W_scat*gevfm;WW=WW2;}

	for(double i=1.;i<=2.;i++)
	{
		for(double j=1.;j<=2.;j++)
		{
			fill_partitionfunction_v2(K[i][j],1*xeps,1.0-1*xeps,anz,phiinv[i][j],phmax[i][j]);
			
			if(ktsplitflag==1)
			{
				fill_partitionfunction(K_u[i][j],u_min,pi,u_anz,phi_u_inv[i][j],phi_u_max[i][j]);
			}
		}
	}

	phix[1.]=phmax[1.][1.];
	phix[2.]=phmax[2.][2.]+phmax[2.][1.];
	probgg=phmax[2.][2.]/phix[2.];

	vector<double> bound=det_boundaries(K[2][2],0.05,1*xeps,1.-xeps);
	b1gg=bound.at(0);
	b2gg=bound.at(1);
	bgg=phmax[2][2];
	make_3c(phiinv[2][2],cgg1,cgg2,cgg3,0.,b1gg,b2gg,bgg);
	bound.clear();

	bqg=phmax[2][1];
	if(bqg>0)
	{
		bound=det_boundaries(K[2][1],0.05,1*xeps,1.-xeps);
		b1qg=bound.at(0);
		b2qg=bound.at(1);
		make_3c(phiinv[2][1],cqg1,cqg2,cqg3,0.,b1qg,b2qg,bqg);
	}
	else
	{
		for(int i=0;i<M;i++){cqg1[i]=0;cqg2[i]=0;cqg3[i]=0;}
	}
	bound.clear();
	bqq=phmax[1][1];
	if(bqq>0)
	{
		bound=det_boundaries(K[1][1],.05,1*xeps,1.-xeps);
		b1qq=bound.at(0);
		b2qq=bound.at(1);
		make_3c(phiinv[1][1],cqq1,cqq2,cqq3,0.,b1qq,b2qq,bqq);
	}
	else
	{
		for(int i=0;i<M;i++){cqq1[i]=0;cqq2[i]=0;cqq3[i]=0;}
	}

	if(ktsplitflag==1)
	{
		make_c(phi_u_inv[2.][1.],cqg_u,0.,phi_u_max[2][1],0*ee_u);
		make_c(phi_u_inv[1.][1.],cqq_u,0.,phi_u_max[1][1],0*ee_u);
		make_c(phi_u_inv[2.][2.],cgg_u,0.,phi_u_max[2][2],0*ee_u);

	}
	srand (time(NULL));
}


double phqq(double z){if(phmax[1][1]>0){return approx_f(z,cqq1,cqq2,cqq3,0.,b1qq,b2qq,bqq);}else{return 0.;}}
double phqg(double z){if(phmax[2][1]>0){return approx_f(z,cqg1,cqg2,cqg3,0.,b1qg,b2qg,bqg);}else{return 0.;}}
double phgg(double z){if(phmax[2][2]>0){return approx_f(z,cgg1,cgg2,cgg3,0.,b1gg,b2gg,bgg);;}else{return 0.;}}
map<double,map<double,fkt>>phiinv1={{1.0, {{1.0,phqq}}},{-1.0, {{-1.0,phqq}}},{2.0, {{1.0,phqg},{2.0,phgg}}}};

double phi_u_invfkt(double z,double i, double j, double (&c)[M])
{
	map<double,double>::iterator it1=phi_u_inv[i][j].begin(),it2=phi_u_inv[i][j].end();
	it2--;
	return approx_f(z, (it1->first),(it2->first),c,ee_u);
}

double phqq_u(double z){return phi_u_invfkt(z,1.0,1.0,cqq_u);}
double phqg_u(double z){return phi_u_invfkt(z,2.0,1.0,cqg_u);}
double phgg_u(double z){return phi_u_invfkt(z,2.0,2.0,cgg_u);}
map<double,map<double,fkt>>phi_u_inv1={{1.0, {{1.0,phqq_u}}},{-1.0, {{-1.0,phqq_u}}},{2.0, {{1.0,phqg_u},{2.0,phgg_u}}}};

double phixg(double x){return phix[2.0]/sqrt(x);}
double phixq(double x){return phix[1.0]/sqrt(x);}
map<double,fkt> phi={{1.0,phixq},{2.0,phixg}};
fvec phi2={phixq,phixg};




double TMDICEevent::tselect(double f,double x, double told)
{
	double colorfact2=1.;
	if(abs(f)==1){colorfact2=colorfact;}
	double r=rand01(),tfin;
	tfin=(told)-log(r)/(phi[abs(f)](x)/(tstar)+colorfact2*W_scat); 
	return tfin;
}

string TMDICEevent::interact_select(double f,double x)
{
	double colorfact2=1.;
	if(abs(f)==1){colorfact2=pow(colorfact,1);}
	double r=rand01();
	if(r<=phi[abs(f)](x)/(phi[abs(f)](x)+tstar*colorfact2*W_scat) ){return "b";}
	if(r>phi[abs(f)](x)/(phi[abs(f)](x)+tstar*colorfact2*W_scat) ){return "c";}
}

double TMDICEevent::typselect(double f1)
{
	double f2;
	double r=rand01();
	if(f1==1.0){f2=1.0;}
	if(f1==-1.0){f2=-1.0;}
	if(f1==2.0)
	{
		if(r<probgg){f2=2.0;}else{double r2=rand01(); if(r2<0.5){f2=1.0;}else{f2=-1.0;} }
	}
	return f2;
}

map<double,map<double,double>> typ2={{1.0,{{1.0,2.0}}},{2.0,{{1.0,-1.0},{-1.0,1.0},{2.0,2.0}} },{-1.0,{{-1.0,2.0}}} };

void  TMDICEevent::fillgen()
{
	vector<TMDICEparticle> newgen,newgenfin;
	if(gen.size()>0)
	{
		for( int j=0; j<gen.size(); j++)
		{
			TMDICEparticle p=gen.at(j);
			
			string sc=interact_select(abs(p.typ),p.x);
			
			if(sc=="b")
			{
				TMDICEparticle p1,p2;
				p1.typ=typselect(p.typ);
				
				p2.typ=typ2[p.typ][p1.typ];
				p1.x=select(phiinv1[abs(p.typ)][abs(p1.typ)],phmax[abs(p.typ)][abs(p1.typ)] ,1*xeps,1.-1*xeps);
				
				p2.x=1.0-p1.x;
				p1.kt=p.kt;
				p2.kt=p1.kt;
				
				p1.kt=p1.kt*p1.x;
				p2.kt=p2.kt*p2.x;

				if(ktsplitflag==1)
				{
					double u_tmp1=0,u_tmp2=0;
					if(phi_u_max[abs(p.typ)][abs(p1.typ)]>0){u_tmp1=select(phi_u_inv1[abs(p.typ)][abs(p1.typ)],phi_u_max[abs(p.typ)][abs(p1.typ)] ,u_min,u_max);}
					double ru1,ru2;

					
					double Q_1=2.0*sqrt(qhatgev*p.x*emax)*sqrt(p1.x*(1.0-p1.x)*f[abs(p.typ)][abs(p1.typ)](p1.x) )*u_tmp1;
					Q_1=sqrt(Q_1);
					double phi_q1=2.0*pi*rand01();
					double Q_2=Q_1;
					double phi_q2=phi_q2+pi;
					if(phi_q2>2.*pi){phi_q2=phi_q2-pi;}
					if(phi_q2<-2.*pi){phi_q2=phi_q2+pi;}
					
					p1.phik=p.phik;
					p2.phik=p.phik;
					
					double ffkt1=1.,ffkt2=1.;
					
					double p1_kx=ffkt1*p1.kt*cos(p1.phik)+Q_1*cos(phi_q1),p1_ky=ffkt1*p1.kt*sin(p1.phik)+Q_1*sin(phi_q1);
					double p2_kx=ffkt2*p2.kt*cos(p2.phik)+Q_2*cos(phi_q2),p2_ky=ffkt2*p2.kt*sin(p2.phik)+Q_2*sin(phi_q2);
					
					
					
					p1.kt=sqrt(p1_kx*p1_kx+p1_ky*p1_ky);
					p1.phik=atan2(p1_ky,p1_kx);

					p2.kt=sqrt(p2_kx*p2_kx+p2_ky*p2_ky);
					p2.phik=atan2(p2_ky,p2_kx);
					
				}
				else
				{
					p1.phik=p.phik;
					p2.phik=p.phik;
				}
				
				p1.x=p1.x*p.x;
				p2.x=p2.x*p.x;
				
				p1.t_old=p.t;
				p2.t_old=p.t;
				p1.t=tselect(abs(p1.typ),p1.x,p.t);
				p2.t=tselect(abs(p2.typ),(p2.x),p.t);
				
				
				if(p1.x>xmin )
				{
					p1.dump=false;
					if(p1.t<tmax){newgen.push_back(p1);}
					if(p1.t>tmax){newgenfin.push_back(p1);}
				}
				else{p1.dump=true;p1.t=tmax;newgenfin.push_back(p1);}
				if(simcasc==true )
				{
					if(p2.x>xmin)
					{
						p2.dump=false;
						if(p2.t<tmax){newgen.push_back(p2);}
						if(p2.t>tmax){newgenfin.push_back(p2);}
					}
					else{p2.dump=true;p2.t=tmax;newgenfin.push_back(p2);}
				}
			}
			
			if(sc=="c")
			{
				double ktnew=WW(rand01());
				ktnew=sqrt(ktnew);
				double phinew=2.*pi*rand01();
				
				TMDICEparticle p1=p;
				double ktold=p1.kt;
				double p1x=p1.kt*cos(p1.phik)+ktnew*cos(phinew),p1y=p1.kt*sin(p1.phik)+ktnew*sin(phinew);
				
				p1.phik=atan2(p1y,p1x);
				p1.kt=sqrt(p1x*p1x+p1y*p1y);

				p1.t_old=p1.t;
				p1.t=tselect(abs(p1.typ),p1.x,p1.t_old);
				
				double p1z=p1.x*emax-(p1.kt*p1.kt)/(4.*p1.x*emax);
				
				if(p1.t<tmax and ktnew>=0*qmin)
				{
					newgen.push_back(p1);
				}
				if( (p1.t>tmax and ktnew>=0*qmin ))
				{
					newgenfin.push_back(p1);
				}
			}
		}
	}
	gen=newgen;
	genfin.insert(genfin.end(),newgenfin.begin(),newgenfin.end());
	casc.insert(casc.end(),gen.begin(),gen.end());
	casc.insert(casc.end(),newgenfin.begin(),newgenfin.end());
}

void TMDICEevent::initgen()
{
	vector<TMDICEparticle>newgen;
	TMDICEparticle p0;
	p0.x=x0;
	p0.typ=typ0;
	p0.t_old=tmin;
	p0.t=tselect(abs(p0.typ),p0.x,tmin);
	p0.kt=kt0;
	p0.phik=0;

	if(p0.x>xmin)
	{
		if(p0.t<tmax){newgen.push_back(p0);}else{genfin.push_back(p0);}
	}
	gen=newgen;
	casc.insert(casc.end(),gen.begin(),gen.end());
	casc.insert(casc.end(),genfin.begin(),genfin.end());
}

void TMDICEevent::make_event()
{
	gen.clear();
	genfin.clear();
	casc.clear();
	initgen();
	while(gen.size()>0)
	{
		fillgen();
	}
}


TMDICEevent::TMDICEevent()
{
	
}

void TMDICEevent::setEmax(double ee)
{
	emax=ee;
}
void TMDICEevent::settmax(double tup)
{
	tmax=tup;
}

void TMDICEevent::settmin(double tdown)
{
	tmin=tdown;
}

void TMDICEevent::setx1(double x00)
{
	x0=x00;
}

void TMDICEevent::setkt1(double kt00)
{
	kt0=kt00;
}

void TMDICEevent::settyp1(double typ00)
{
	typ0=typ00;
}

TMDICEevent:: ~TMDICEevent()
{
	
}


