#! /bin/bash

rm -rf parmetis-4.0.3
wget http://ppmcore.mpi-cbg.de/upload/parmetis-4.0.3.tar.gz -O parmetis-4.0.3.tar.gz
tar -xf parmetis-4.0.3.tar.gz
cd parmetis-4.0.3
# Change to 64 bit parmetis

if [[ "$OSTYPE" == "darwin"* ]]; then
  sed_command=gsed
else
  sed_command=sed
fi

$sed_command -i "/#define\sIDXTYPEWIDTH\s32/c\#define IDXTYPEWIDTH 64" metis/include/metis.h

CFLAGS=-fPIC make config prefix=$1/PARMETIS
make -j $2
if [ $? -ne 0 ]; then
    echo "PARMETIS error installing"
    exit 0
fi
mkdir $1/PARMETIS
make install

#### Apply patch if we are on cygwin

if [[ "$OSTYPE" == "cygwin" ]]; then
  cd $1/PARMETIS/include
  wget http://openfpm.mpi-cbg.de/upload/parmetis_patch -O parmetis_patch
  patch < parmetis_patch
fi


