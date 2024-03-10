#include<fstream>
#include<string>
#include<iostream>
#include<iomanip>
#include "vle.h"
#include "TMDICE.h"
#include "functions_lib.h"

using namespace std;

int main(int argc, char **argv)
{

	string inputfile="inputfile";
	readTMDICEparameters(inputfile);
	double Nevents=1*pow(10,7);
	setTMDICE();
	ofstream o;
	o.open(argv[1]);
	o<<"i,x,typ"<<endl;

	ifstream s;
	s.open(inputfile.c_str());
	while(s.eof()==false)
	{
		string s_in;
		getline(s,s_in);
	}
	s.close();

	for(int i=0;i<Nevents;i++)
	{
		TMDICEevent jet;
		jet.cascfill=false;
		jet.filltemps=false;
		cout<<static_cast<int>(100*static_cast<double>(i)/Nevents)<<"%\r";
		jet.dumpthedump=true;
		jet.setx1(1.0);
		
		jet.make_event();;
		vector<TMDICEparticle*> partlist=jet.getgenfin();

		for(int j=0; j<partlist.size();j++)
		{
			TMDICEparticle *part=partlist.at(j);			
			if(part->dump==false ){o<<i<<","<<part->x<<","<<part->typ<<endl;}
		}
		partlist={};
		jet.genfin={};
	}
	o.close();
	cout<<endl<<"Program finished."<<endl;
	cout<<"Output written to file: "<<argv[1]<<endl;
	
}
