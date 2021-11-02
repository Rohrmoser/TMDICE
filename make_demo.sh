#!/bin/bash

if [ $# -ne 1 ]
  then
    echo "Please type a path after 'bash make_demo.sh' "
fi

mkdir -p $1
cp demo.cpp $1/
cp inputfile $1/
cp evaluation_demo.nb $1/

echo "#!/bin/bash" > $1/"makefile.sh"
echo "g++ -std=c++11 ""demo.cpp "$PWD"/TMDICE.cpp "$PWD"/TMDICE_lib.cpp -I"$PWD" -O3 -o ""demo.out">>$1/"makefile.sh"
echo "echo \"Created the TMDICE demo. To run it, replace in the following line <outputfile> with a filename of your choice and type the resulting line in a terminal window, then press enter: \" "'$PWD'"\"/demo.out <outputfile>\"">>$1/"makefile.sh"

echo "Created the demo  in " $1
echo "To compile the demo  type: cd "$1" ; bash makefile.sh" 
