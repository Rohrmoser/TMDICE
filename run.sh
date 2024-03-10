#!/bin/bash
cc="g++"
inputdir=./
$cc runtmdice.cpp $inputdir/globalconsts.cpp $inputdir/TMDICE.cpp $inputdir/BDIM.cpp $inputdir/vle.cpp $inputdir/TMDICE_base.cpp $inputdir/TMDICE_lib.cpp $inputdir/decor.cpp $inputdir/functions_lib.cpp $inputdir/medkernels.cpp  -I$inputdir/ -O3 -o a.out

