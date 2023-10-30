#!/bin/bash

# Set the base directory to the directory where the script is located
cd "$(dirname "$0")"

# Create build folder
mkdir build && cd build 

# Check if an argument was provided
if [ -n "$1" ]; then
  if [ -n "$2" ]; then
    if [ -n "$3" ]; then
	  # Run cmake with the provided path as the MUI include and lib directories
	  cmake -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=../ .. -DMUI_C_INCLUDE_DIR=$1 -DMUI_C_LIB_DIR=$2 -DCMAKE_PREFIX_PATH=$3
    else
	  # Run cmake with the provided path as the MUI include and lib directories
	  cmake -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=../ .. -DMUI_C_INCLUDE_DIR=$1 -DMUI_C_LIB_DIR=$2
    fi
  else
    # Run cmake with the provided path as the MUI base directory
    cmake -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=../ .. -DMUI_BASE_DIR=$1 -DCMAKE_PREFIX_PATH=$1/share/MUI-2.0/cmake
  fi
else
  # Run cmake with the default MUI directory
  cmake -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=../ ..
fi

# Run make to build the executable
make 2>&1 | tee make.log && cd .. 

START=$SECONDS

mpirun -np 1 heat-right :\
       -np 1 heat-left | tee output.log

DURATION=$(( SECONDS - START ))

if (( $DURATION > 3600 )) ; then
    let "hours=DURATION/3600"
    let "minutes=(DURATION%3600)/60"
    let "seconds=(DURATION%3600)%60"
    echo "Completed in $hours hour(s), $minutes minute(s) and $seconds second(s)"
elif (( $DURATION > 60 )) ; then
    let "minutes=(DURATION%3600)/60"
    let "seconds=(DURATION%3600)%60"
    echo "Completed in $minutes minute(s) and $seconds second(s)"
else
    echo "Completed in $DURATION seconds"
fi

paraview Resources/view_per_iteration.pvsm