#! /bin/bash

set -euo pipefail

rm -rf ParMETIS
wget http://ppmcore.mpi-cbg.de/upload/ParMETIS_master.tar.gz -O ParMETIS_master.tar.gz
tar -xf ParMETIS_master.tar.gz
cd ParMETIS
make distclean
# perl -ni.old -e 'print;if ($.==9) {print "target_link_libraries(metis -lGKlib)\n"}' libmetis/CMakeLists.txt
CMAKE_POLICY_VERSION_MINIMUM=3.5 make config shared=1 cc=$3 prefix=$1/PARMETIS gklib_path=$1/METIS metis_path=$1/METIS
make -j $2
make install
rm -f "$1/PARMETIS/lib/libparmetis.a"
