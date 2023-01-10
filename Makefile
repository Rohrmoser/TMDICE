inputdir=../src/


cc=g++ -std=c++11
vle:  
	g++ -std=c++11 $(inputdir)demovle.cpp $(inputdir)globalconsts.cpp $(inputdir)vle.cpp $(inputdir)TMDICE_base.cpp $(inputdir)TMDICE_lib.cpp $(inputdir)decor.cpp  -I$(inputdir) -O3 -o demovle.out

bdim: 
	$(cc) $(inputdir)demobdim.cpp $(inputdir)globalconsts.cpp $(inputdir)BDIM.cpp $(inputdir)TMDICE_base.cpp $(inputdir)TMDICE_lib.cpp $(inputdir)decor.cpp $(inputdir)functions_lib.cpp $(inputdir)medkernels.cpp -I$(inputdir) -O3 -o demobdim.out

tmdice: 
	$(cc) $(inputdir)demo.cpp $(inputdir)globalconsts.cpp $(inputdir)TMDICE.cpp $(inputdir)BDIM.cpp $(inputdir)vle.cpp $(inputdir)TMDICE_base.cpp $(inputdir)TMDICE_lib.cpp $(inputdir)decor.cpp $(inputdir)functions_lib.cpp $(inputdir)medkernels.cpp -I$(inputdir) -O3 -o demo.out
	



