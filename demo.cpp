#include<fstream>
#include<string>
#include<iostream>
#include "TMDICE.h"

int main(int argc, char **argv)
{
	string inputfile="inputfile";
	readTMDICEparameters(inputfile);

	double Nevents=pow(10,5);
	setTMDICE();

	
	ofstream o;
	o.open(argv[1]);
	
	o<<"Parameters"<<endl;
	ifstream s;
	s.open(inputfile.c_str());
	while(s.eof()==false)
	{
		string s_in;
		getline(s,s_in);
		o<<s_in<<endl;
	}
	s.close();
	o<<endl<<"kt...momentum component transverse to jet-axis"<<endl;
	o<<"phik...azimuthal angle of transverse momentum component"<<endl;
	o<<"typ...1 for quarks, 2 for gluons"<<endl;
	o<<endl;
	o<<"Index for event"<<"	"<<"x"<<"	"<<"kt"<<"	"<<"phik"<<"	"<<"typ"<<endl;
	for(int i=0;i<Nevents;i++)
	{
		if(100*i/Nevents==round(100*i/Nevents)){cout<<100*i/Nevents<<"% of events generated"<<endl;}
		TMDICEevent jet;
		jet.make_event();
		
		for(int j=0;j<jet.genfin.size();j++)
		{
			if(jet.genfin.at(j).dump==false){o<<i<<" "<<jet.genfin.at(j).x<<" "<<jet.genfin.at(j).kt<<" "<<jet.genfin.at(j).phik<<" "<<jet.genfin.at(j).typ<<endl;}
		}
		o<<"new event"<<endl;
	}
	o.close();
	cout<<"Program finished."<<endl;
	cout<<"Output written to file: "<<argv[1]<<endl;
	
}