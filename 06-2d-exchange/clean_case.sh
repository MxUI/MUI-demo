#!/bin/bash

# Set the base directory to the directory where the script is located
cd "$(dirname "$0")"

# Clean executables
cd build && make clean && cd ..

# Clean log and results
rm -f make.log
rm -f output.log
rm -f core.*
rm -rf build
rm -f solution_coarse002.dat
rm -f solution_coarse004.dat
rm -f solution_refine002.dat
rm -f solution_refine004.dat

