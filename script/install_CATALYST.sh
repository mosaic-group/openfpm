#! /bin/bash

rm -rf catalyst
git clone https://gitlab.kitware.com/paraview/catalyst.git
cd catalyst
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$1/CATALYST ..
make -j $2
make install