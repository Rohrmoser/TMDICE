#include<fstream>
#include<string>
#include<iostream>
#include<iomanip>
#include "vle.h"


int main(int argc, char **argv)
{
	string inputfile="inputfilevle";
	readTMDICEparameters(inputfile);

	double Nevents=1*pow(10,6);
	setVLE();

	ofstream o;
	o.open(argv[1]);
	//o<<"ind,x,xold,q,qold,kt12,th12,th12old,adr"<<endl;
	//o<<"ind,typ,x,kt,th12,pt,eta,phi,adr"<<endl;
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
	o<<endl;
	o<<"Index for event"<<"	"<<"x"<<"	"<<"kt"<<"	"<<"phik"<<"	"<<"typ"<<endl;*/
	loadstart=false;
	for(int i=0;i<Nevents;i++)
	{
//		cout<<i<<endl;
		//if(100*i/Nevents==round(100*i/Nevents)){cout<<"\r\r"/*"\033[F"*/<<100*i/Nevents<<"% of events generated"<<endl;}
		loadingbar(0,static_cast<double>(Nevents),static_cast<double>(i),Nevents);
		VLEevent jet(stop_non_resolved_but_resolvable,stop_non_resolved_but_resolvable);//(stop_when_td_exceeds_tmax_and_ktmin,stop_when_td_exceeds_tmax_and_ktmin,bias_tdtl);
		
//double x0=select_gausscdf();
//		cout<<"x0="<<x0<<endl;
		jet.setx1(1.0);
		
		jet.make_event();//2VLEinit();
		//jet.initgen();
		//jet.fillgen();
		//jet.fillgen();
		

		for(int j=0;j<jet.genfin.size();j++)
		{
			if(jet.genfin.at(j).dump==false){o<<i<<" "<<jet.genfin.at(j).x<<" "<<jet.genfin.at(j).kt<<" "<<jet.genfin.at(j).phik<<" "<<jet.genfin.at(j).typ<<endl;}
		}
		
	//	o<<"new event"<<endl;
	}
	o.close();
	cout<<endl<<"Program finished."<<endl;
	cout<<"Output written to file: "<<argv[1]<<endl;
	
}