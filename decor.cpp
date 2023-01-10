#include "deps.h"
#include "decor.h"

bool func_exec_tmdice=false;
bool loadstart=false;


void disclaimer()
{
	if(func_exec_tmdice==false)
	{
	cout<<endl;
	cout<<"This is the \033[1;31mTMDICE\033[0m program"<<endl;
	cout<<"by Martin Rohrmoser, 2021."<<endl;
	cout<<"If you use it, please cite it!"<<endl<<endl;
	}
	func_exec_tmdice=true;
}

void loadingbar(double valmin,double valmax,double val,double valanz)
{
	if(loadstart==false)
	{
		cout<<"1%        10%       20%       30%       40%       50%       60%       70%       80%       90%      100%"<<std::flush<<endl;
		cout<<"|        |         |         |         |         |         |         |         |         |         |"<<std::flush<<endl;
		//     [****************************************************************************************************]
		//      ****************************************************************************************************
		loadstart=true;
	}
	//double prog=100*(val-valmin)/(valmax-valmin);

	//double prog=static_cast<int>(floor(100*(val-valmin)/(valmax-valmin)))-static_cast<int>(floor(100*(val-dt-valmin)/(valmax-valmin)));
	//if(prog>0)
	//if(prog==round(prog))
	double procent=(valmax-valmin)/100;
	double ww=(val-valmin)/procent;
//	if(ww-floor(ww)<=0.01/procent)
if(ww>=floor(ww) and ww-(valmax-valmin)/(valanz*procent)<floor(ww))
//if(ww==floor(ww))
	{
		cout<<"*"<<std::flush;
	}
	/*if(abs(ww-100)<=0.00001)
	{
		cout<<std::flush<<endl;
		loadstart=false;
	}*/
}
