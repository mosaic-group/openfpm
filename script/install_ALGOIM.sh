#! /bin/bash

wget http://ppmcore.mpi-cbg.de/upload/algoim.tar.gz -O algoim.tar.gz
tar -xf algoim.tar.gz
mv algoim $1/ALGOIM
mv $1/ALGOIM/src $1/ALGOIM/include
