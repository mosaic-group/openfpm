#! /bin/bash

rm -rf minter
wget http://ppmcore.mpi-cbg.de/upload/minter.tar.gz -O minter.tar.gz
tar -xf minter.tar.gz

mkdir -p $1/MINTER
mv minter/include $1/MINTER/minter
