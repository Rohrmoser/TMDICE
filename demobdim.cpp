#include<fstream>
#include<string>
#include<iostream>
#include "BDIM.h"

map<double,double> gausscdf;

void fill_gausscdf(double sig)
{
	double anz=10000.;
	double dx=1./anz;
	double norm=0,normtmp=0;
	for(double x=0;x<=1;x=x+dx)
	{
		norm=norm+dx*exp(-pow(1.-x,2)/(2.*sig*sig));
	}
	
	for(double x=0;x<=1;x=x+dx)
	{
		normtmp=normtmp+dx*exp(-pow(1.-x,2)/(2.*sig*sig));
		gausscdf[normtmp/norm]=x+dx;
	}
}

double select_gausscdf()
{
	double r=rand01();
	map<double,double>::iterator it=gausscdf.lower_bound(r);
	return it->second;
}

int main(int argc, char **argv)
{
	string inputfile="inputfilebdim";
	readTMDICEparameters(inputfile);

	double Nevents=pow(10,6);
	setBDIM();
fill_gausscdf(pow(10,-2));
	
	ofstream o;
	o.open(argv[1]);
	
//	o<<"Parameters"<<endl;
	ifstream s;
	s.open(inputfile.c_str());
	while(s.eof()==false)
	{
		string s_in;
		getline(s,s_in);
//		o<<s_in<<endl;
	}
	s.close();
/*	o<<endl<<"kt...momentum component transverse to jet-axis"<<endl;
	o<<"phik...azimuthal angle of transverse momentum component"<<endl;
	o<<"typ...1 for quarks, 2 for gluons"<<endl;
	o<<endl;*/
	o<<"Index x kt phik typ"<<endl;
	loadstart=false;
	for(int i=0;i<Nevents;i++)
	{
//		cout<<i<<endl;
		//if(100*i/Nevents==round(100*i/Nevents)){cout<<"\r\r"/*"\033[F"*/<<100*i/Nevents<<"% of events generated"<<endl;}
		loadingbar(0,static_cast<double>(Nevents),static_cast<double>(i),Nevents);
		BDIMevent jet;
		
		//double x0=select_gausscdf();
//		cout<<"x0="<<x0<<endl;
		jet.setx1(1.0);
		
		jet.make_event();
		
		for(int j=0;j<jet.genfin.size();j++)
		{
			if(jet.genfin.at(j).dump==false /*and jet.genfin.at(j).islead==1*/)
			{o<<i<<" "<<jet.genfin.at(j).x<<" "<<jet.genfin.at(j).kt<<" "<<jet.genfin.at(j).phik<<" "<<jet.genfin.at(j).typ<<endl;}
		}
	//	o<<"new event"<<endl;
	}
	o.close();
	cout<<endl<<"Program finished."<<endl;
	cout<<"Output written to file: "<<argv[1]<<endl;
	
}