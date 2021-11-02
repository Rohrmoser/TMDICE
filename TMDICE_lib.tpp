#include<iostream>
#include<math.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "TMDICE_lib.h"
using namespace std;

	

   	template<typename T, typename ...Rest> double integral(  double f(  T x,Rest...rest), double a,  double b,T x, Rest...rest )
    	{
    	double iwert=0;
		if(b>a){
			for(int i=1;i<=n;i++)
			{
			iwert=iwert+
			f((a+b)/2.0f+((b-a)/2.0f)*gaussx[i], rest...)*gaussw[i];
			}
		}
	return iwert*pi/n*(b-a)/2.0f;
	}
