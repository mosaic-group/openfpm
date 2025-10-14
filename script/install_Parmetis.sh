#! /bin/bash

rm -rf ParMETIS
wget http://ppmcore.mpi-cbg.de/upload/ParMETIS_master.tar.gz -O ParMETIS_master.tar.gz
tar -xf ParMETIS_master.tar.gz
cd ParMETIS
make distclean
# perl -ni.old -e 'print;if ($.==9) {print "target_link_libraries(metis -lGKlib)\n"}' libmetis/CMakeLists.txt
make config shared=1 cc=$3 prefix=$1/ParMETIS gklib_path=$1/METIS metis_path=$1/METIS
make -j $2
make install
