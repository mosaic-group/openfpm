#! /bin/bash

prefix_dependencies="$1"
prefix_openfpm="$2"

if [[ "$OSTYPE" == "linux-gnu" || "$OSTYPE" == "linux" || "$OSTYPE" == "cygwin" ]]; then
  bash_library="export LD_LIBRARY_PATH=\"\$LD_LIBRARY_PATH:/$prefix_openfpm/openfpm_devices/lib:/$prefix_openfpm/openfpm_vcluster/lib"
else
  bash_library="export DYLD_LIBRARY_PATH=\"\$DYLD_LIBRARY_PATH:/$prefix_openfpm/openfpm_devices/lib:/$prefix_openfpm/openfpm_vcluster/lib"
fi

bash_path="export PATH=\""

if [ -d "$prefix_dependencies/MPI" ]; then
  bash_path="$bash_path$prefix_dependencies/MPI/bin:"
  bash_library="$bash_library:$prefix_dependencies/MPI/lib"
fi

if [ -d "$prefix_dependencies/METIS" ]; then
  bash_library="$bash_library:$prefix_dependencies/METIS/lib"
fi

if [ -d "$prefix_dependencies/PARMETIS" ]; then
  bash_library="$bash_library:$prefix_dependencies/PARMETIS/lib"
fi

if [ -d "$prefix_dependencies/BOOST" ]; then
  bash_library="$bash_library:$prefix_dependencies/BOOST/lib"
fi

if [ -d "$prefix_dependencies/HDF5/lib" ]; then
  bash_library="$bash_library:$prefix_dependencies/HDF5/lib"
elif [ -d "$prefix_dependencies/HDF5/lib64" ]; then
  bash_library="$bash_library:$prefix_dependencies/HDF5/lib64"
fi

if [ -d "$prefix_dependencies/LIBHILBERT" ]; then
  bash_library="$bash_library:$prefix_dependencies/LIBHILBERT/lib"
fi

if [ -d "$prefix_dependencies/PETSC" -a -f "$prefix_dependencies/PETSC/include/petsc.h" ]; then
  bash_library="$bash_library:$prefix_dependencies/PETSC/lib"
fi

if [ -d "$prefix_dependencies/OPENBLAS" ]; then
  bash_library="$bash_library:$prefix_dependencies/OPENBLAS/lib"
fi

if [ -d "$prefix_dependencies/SUITESPARSE"  -a -f "$prefix_dependencies/SUITESPARSE/include/umfpack.h" ]; then
  bash_library="$bash_library:$prefix_dependencies/SUITESPARSE/lib"
fi

if [ -d "$prefix_dependencies/CATALYST/lib" ]; then
  bash_library="$bash_library:$prefix_dependencies/CATALYST/lib"
elif [ -d "$prefix_dependencies/CATALYST/lib64" ]; then
  bash_library="$bash_library:$prefix_dependencies/CATALYST/lib64"
fi

# in cygwin we have to add to PATH additional directories
if [[ "$OSTYPE" == "cygwin" ]]; then
	bash_path="$bash_path:$prefix_dependencies/BOOST/bin:$prefix_dependencies/HDF5/bin"
fi

bash_path="$bash_path:\$PATH\""
bash_library="$bash_library\""

##### Writing openfpm_vars file

echo "$bash_path" > openfpm_vars
echo "$bash_library" >> openfpm_vars
echo "export PURE_PYTHON=1" >> openfpm_vars