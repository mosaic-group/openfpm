#! /bin/bash

rm -rf GKlib
wget http://ppmcore.mpi-cbg.de/upload/GKlib_master.tar.gz -O GKlib_master.tar.gz
tar -xf GKlib_master.tar.gz
cd GKlib
make distclean
make config cc=$3 shared=1 prefix=$1/METIS
make -j $2
make install
cd ..

rm -rf METIS
wget http://ppmcore.mpi-cbg.de/upload/METIS_master.tar.gz -O METIS_master.tar.gz
tar -xf METIS_master.tar.gz
cd METIS
make distclean
perl -ni.old -e 'print;if ($.==9) {print "target_link_libraries(metis -lGKlib)\n"}' libmetis/CMakeLists.txt
if [ ! -f $1/METIS/lib/libGKlib.so  ]; then
	if [ -f $1/METIS/lib/libGKlib.so.0  ]; then
		ln -s $1/METIS/lib/libGKlib.so.0  $1/METIS/lib/libGKlib.so
	fi
fi

make config shared=1 i64=1 r64=1 cc=$3 prefix=$1/METIS gklib_path=$1/METIS
make -j $2
make install
