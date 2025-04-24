#! /bin/bash

rm -rf minter
git clone -b header_only https://git.mpi-cbg.de/mosaic/software/math/minter.git

mkdir -p $1/MINTER
mv minter/include $1/MINTER/minter