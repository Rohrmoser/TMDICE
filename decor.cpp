#include "deps.h"
#include "decor.h"
using namespace std;


bool func_exec_tmdice=false;

void disclaimer()
{
	if(func_exec_tmdice==false)
	{
	cout<<endl;
	cout<<"This is the \033[1;31mTMDICE\033[0m program\n";
	cout<<"by Martin Rohrmoser\n";
	cout<<"If you use it, please cite it!\n\n";
	}
	func_exec_tmdice=true;
}

