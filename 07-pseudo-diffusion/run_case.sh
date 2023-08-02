#!/bin/bash

# Set the base directory to the directory where the script is located
cd "$(dirname "$0")"

# Create build folder
mkdir build && cd build 

# Check if an argument was provided
if [ -n "$1" ]; then
  # Run cmake with the provided path as the MUI directory
  cmake -DCMAKE_PREFIX_PATH=$1 -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=../ ..
else
  # Run cmake with the default MUI directory
  cmake -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=../ ..
fi

# Run make to build the executable
make 2>&1 | tee make.log && cd .. 

START=$SECONDS

mpirun -np 2 3D-pseudo-diffusion-coarse : \
       -np 2 3D-pseudo-diffusion-fine 2>&1 | tee output.log

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

#paraview Resources/view.pvsm
