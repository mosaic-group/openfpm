#!/bin/bash 

rm -rf boost_1_84_0
wget http://ppmcore.mpi-cbg.de/upload/boost_1_84_0.tar.gz -O boost_1_84_0.tar.gz
tar -xf boost_1_84_0.tar.gz
cd boost_1_84_0

./bootstrap.sh --with-toolset=$3

TOOLSET=""

if [ x"$CXX" == x"icpc" ]; then
	TOOLSET="intel-linux"
elif [ x"$CXX" == x"clang++" ]; then
	TOOLSET="clang"
else
	if [[ x"$CXX" == x"g++"* ]]; then
		g++ --version | grep "Apple LLVM\|Apple clang" >/dev/null 2>&1
		if [ $? == 0 ]; then
			TOOLSET="clang"
		fi
	else
		TOOLSET="gcc"
	fi
fi

mkdir $1/BOOST
# Several flavours
if [[ x"$OSTYPE" == x"darwin"* ]]; then
    if [ x$(uname -m) == x"arm64" ]; then
	./b2 -a -j $2 install --prefix=$1/BOOST address-model=64 architecture=arm abi=aapcs binary-format=mach-o toolset=$TOOLSET  -sNO_LZMA=1 -sNO_ZSTD=1
    else
	./b2 -a -j $2 install --prefix=$1/BOOST address-model=64 architecture=x86 abi=sysv binary-format=mach-o toolset=$TOOLSET  -sNO_LZMA=1 -sNO_ZSTD=1
    fi
else
    ./b2 -a -j $2 install --prefix=$1/BOOST  -sNO_LZMA=1 -sNO_ZSTD=1
fi
